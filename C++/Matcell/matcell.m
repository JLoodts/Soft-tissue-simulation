% EEM_Staaf_1
%
% Eindige Elementen Methode toegepast op 2 dimensionale staafmodellen  
%
% script geschreven voor de oefening HJ35:
%    Numerieke modellering in de mechanica
%
clear all; close all;

% ----------------------------------------------------------------------------- %
% INITIALISATIE
% ----------------------------------------------------------------------------- %

Node = [%   nr     xpos    ypos
            1,     0,      0;
            2,     2.5,    4.33013;
            3,     7.5,    4.33013;
            4,     10,     0;
            5,     5,      0        ];
BeginSurface = getSurface(Node)
     
Element = [% nr     begin   end
             1,     1,      2;
             2,     2,      3;
             3,     3,      4;
             4,     1,      5;
             5,     2,      5;
             6,     3,      5;
             7,     4,      5       ];
     
% initialisatie van stijfheidsmatrix, rechterlid en oplossing
%    nElem = aantal elementen
%    nNode = aantal knooppunten
%    K     = totale stijfheidsmatrix
%    F     = rechterlid
%    d     = oplossing
 nElem = 7; 
 nNode = 5;
 nDOF  = 2*nNode;
 K     = zeros(nDOF,nDOF);
 f     = zeros(nDOF,1);
 d     = zeros(nDOF,1);

% ----------------------------------------------------------------------------- %
% OPBOUW VAN TOTALE STIJFHEIDSMATRIX
% ----------------------------------------------------------------------------- %

% lus over elementen
 for p = 1:nElem

   A    = 0.003;                % de dwarsdoorsnede [m²]
   E    = 200E3;                % elasticiteitsmodulus [N/m²]
   pid1 = Element(p,2);         % ID van het beginpunt van het p-de element
   pid2 = Element(p,3);         % ID van het eindpunt van het p-de element
   x1   = Node(pid1,2);         % x-coordinaat van het beginpunt
   y1   = Node(pid1,3);         % y-coordinaat van het beginpunt
   x2   = Node(pid2,2);         % x-coordinaat van het eindpunt
   y2   = Node(pid2,3);         % y-coordinaat van het eindpunt
   L    = sqrt((x2-x1)^2+(y2-y1)^2); % initiele lengte van de p-de veer
   
   % stel elementstijfheidsmatrix Ke_loc op in lokaal assenstelsel
   Ke_loc = -(E*A/L)*...
            [ -1  0  1  0;
               0  0  0  0;
               1  0 -1  0;
               0  0  0  0  ];
    
   % stel transformatiematrix op
   cosfi = (x2-x1)/L;
   sinfi = (y2-y1)/L;
   G = [ cosfi -sinfi     0      0;
         sinfi  cosfi     0      0;
             0      0 cosfi -sinfi;
             0      0 sinfi  cosfi  ];
   
   % stel elementstijfheidsmatrix op in globaal assenstelsel
   Ke_glo = G*Ke_loc*G';
   
   % assemblage
   pid1_glo = pid1*2-1;
   pid2_glo = pid2*2-1;
   index = [pid1_glo, pid1_glo+1, pid2_glo, pid2_glo+1];
   K(index, index) = K(index, index) + Ke_glo;
    
end % van loop over elementen

% visualisatie van de totale stijfheidsmatrix
% spy(K);
 
 % ----------------------------------------------------------------------------- %
% OPBOUW VAN RECHTERLID
% ----------------------------------------------------------------------------- %

% applying BC: known forces
% op node 4 werkt een kracht in de x-richting [N]
%f(3*2-1) =  100000;
 
 % ----------------------------------------------------------------------------- %
% PARTITIONERING
% ----------------------------------------------------------------------------- %

% index_c = vrijheidsgraden in d, die beperkt zijn
% index_f = vrijheidsgraden in d, die onbekend zijn
% K_ff    = deel van stijheidsmatrix geassocieerd met onbekende vrijheidsgraden
% K_cf    = deel van stijfheidsmatrix nodig voor berekening van reactiekrachten
% F_f     = deel van rechterlid met bekende belasting
index_c = [ 1*2-1,      % node 1 x-direction   
            1*2,        % node 1 y-direction
            4*2-1,      % node 4 x-direction
            4*2,        % node 4 y-direction
                    ];   
d(index_c) = [ -1,      % node 1 x-direction   
               0,        % node 1 y-direction
               1,      % node 4 x-direction
               0,        % node 4 y-direction
                    ];  
index_f = setdiff( 1:nDOF, index_c );
K_ff    = K(index_f, index_f);
K_cf    = K(index_c, index_f);
f_bc    = f;
% applying BC: known displacements as defined in d(index_c)
for i = 1:length(index_c),
    f_bc = f_bc - d(index_c(i))*K(:,index_c(i)); 
end
f_f     = f_bc(index_f);

% ----------------------------------------------------------------------------- %
% OPLOSSING
% ----------------------------------------------------------------------------- %
 
% oplossen van het stelsel vergelijkingen naar de onbekende vrijheidsgraden
d_f = K_ff\f_f;
 
% berekenen van reactiekrachten
f_c = K_cf*d_f;

% zet de berekende vrijheidsgraden op hun plaats in d
d(index_f) = d_f;

% sla oplossing op in Displacement matrix
 Displacement = zeros(size(Node));
 Displacement(:,1) = Node(:,1);
 Displacement(:,2) = d(1:2:10,1);
 Displacement(:,3) = d(2:2:10,1);

 DisplacedNode = Node + Displacement;
 EndSurface = getSurface(DisplacedNode)
 
% laat de resultaten zien
figure;
plot(Node(:,2), Node(:,3),'bo');
hold on;
plot(DisplacedNode(:,2), DisplacedNode(:,3),'r*');
axis equal;
legend('original configuration','configuration under load');
% exporteer resultaat naar Femap neutral file
%  NeuName = 'vakwerk7.neu';
%  NeuWrite( NeuName, Displacement, [], [] );

% ----------------------------------------------------------------------------- %
% POSTPROCESSING
% ----------------------------------------------------------------------------- %

% initialisatie
 RodForce = zeros( nElem, 3 ); % staafkrachtenresultaat voor NeuWrite
 
% lus over elementen
 for p = 1:nElem
    
 end
 
% exporteer resultaat naar Femap neutral file
%  NeuName = 'vakwerk7.neu';
%  NeuWrite( NeuName, Displacement, RodForce, [] );