function y = periospline(x,f,t)

% x: vector van n+1 abscissen
% f: [d x n] matrix die de functiewaarden van d verschillende functies bevat
% t: de N evaluatiepunten
% y: [d x N] matrix met in elke rij de functiewarden in de punten t voor de overeenkomstige rij van f

ttime = clock;

[d,n] = size(f);        % n = aantal interpolatiepunten - 1
fmod = [f,f(:,1)];      % om straks makkelijker df te berekenen
N = length(t);          % N = aantal avaluatiepunten

%--------------------------------------------------------------------------
% bepaling van sdd: de tweede afgeleide van de spline s in de punten xi
%--------------------------------------------------------------------------

% enkele differenties
dx = x(2:n+1)-x(1:n); dx = [dx,dx(1)];
df = fmod(:,2:n+1)-fmod(:,1:n);

% berekenen van de coeffientenmatrix A
A = zeros(n,n);
A(1,1) = 2*(dx(n)+dx(1)); A(1,2) = dx(1); A(1,n) = dx(n);
for(i=2:n-1)
    A(i,i-1) = dx(i-1);
    A(i,i)   = 2*(dx(i-1)+dx(i));
    A(i,i+1) = dx(i);
end
A(n,1)   = dx(n); A(n,n-1) = dx(n-1); A(n,n)   = 2*(dx(n-1)+dx(n));

% het rechterlid F
F(:,1) = df(:,1)/dx(1)-df(:,n)/dx(n);
for(i=2:n)
    F(:,i) = df(:,i)/dx(i)-df(:,i-1)/dx(i-1);
end
F = 6*F';

% oplossen van het stelsel
sdd = A\F;
sdd = [sdd;sdd(1,:)]; sdd = sdd';

%sdd = 2*ones(d,n+1);

%--------------------------------------------------------------------------
% berekening van de functiewaarden y = s(t(i)) aan de hand van sdd
%--------------------------------------------------------------------------

% herschalen van t zodat alle t(i) binnen [xmin,xmax] vallen
span = x(n+1)-x(1); tmod = mod(t-x(1),span); tmod = tmod + x(1);

% enkele tellers voor de while lussen
i = 2; j = 1; nrreset = 0; % nrreset om bij foutieve invoer oneindige lus te vermijden
while ((j < N)&(nrreset<N))
    while (~((x(i-1)<= tmod(j))&(tmod(j) < x(i))))&(nrreset<N)
        i = i+1;
        if(i>=n+2) i = 2; nrreset = nrreset+1; end
    end
    xm_xim1 = tmod(j)-x(i-1);
    xim_x   = x(i)-tmod(j);
    y(:,j) = (fmod(:,i)*xm_xim1+fmod(:,i-1)*xim_x)/dx(i) ...
        + 1/6*(xm_xim1^3/dx(i)-dx(i)*xm_xim1)*sdd(:,i) ...
        + -1/6*(-xim_x^3/dx(i)+dx(i)*xim_x)*sdd(:,i-1) ;
    j = j+1;
end
etime(clock,ttime)