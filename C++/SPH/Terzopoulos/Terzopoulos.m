% Terzopoulos Demetri, Platt John, Barr alan, Fleischer Kurt
% Elastically Deformable Models, 1987,
% Computer Graphics, Vol 21:4, p. 205-214.

function [Movie] = Terozpoulos()

close all;
clear all;
% r[a1,a2] maps the points from parameter space [a1,a2] to the real world space (r1,r2)
% 0 <= a1,a2 <= 1
constM = 25; % nuconstMber of grid points in the a1-direction
constN = 20; % constNumber of grid points in the a2-direction

x = linspace(0,1,constN);
y = linspace(0,1,constM);
h1 = y(2)-y(1); % inter-node spacing in the a1-direction
h2 = x(2)-x(1); % inter-node spacing in the a2-direction
[xx,yy] = meshgrid(x,y); 
R0x = MatrixToVector(xx);
R0y = MatrixToVector(yy);
clear x y xx yy;

Rt_x = R0x; % the result of the mapping
Rt_y = R0y; % the result of the mapping
Rtminusdt_x = Rt_x;
Rtminusdt_y = Rt_y;

Rx = VectorToMatrix(Rt_x,constM,constN);
Ry = VectorToMatrix(Rt_y,constM,constN);  
Kt = CalculateK(Rx,Ry,constM,constN,h1,h2);
M = CalculateM(constM,constN);
C = CalculateC(constM,constN);
[F_x,F_y] = CalculateF(constM,constN);


% integration constants
tmin = 0;
tmax = 1;
dt = 0.1;

MandC = (1/dt^2)*M + (1/(2*dt))*C;

% movie 
XMIN = -2.1; XMAX = 2.1;
YMIN = -4; YMAX = 1;
frameNr = 1;
mapToSpace(Rx,Ry);
AXIS([XMIN XMAX YMIN YMAX]);
Movie(frameNr) = getframe;

for t = tmin:dt:tmax
    
Rx = VectorToMatrix(Rt_x,constM,constN);
Ry = VectorToMatrix(Rt_y,constM,constN);    
Kt = CalculateK(Rx,Ry,constM,constN,h1,h2);

Vt_x = (Rt_x - Rtminusdt_x)/dt;
Vt_y = (Rt_y - Rtminusdt_y)/dt;

At = Kt + MandC;
Gt_x = F_x + MandC*Rt_x + ((1/dt)*M - (1/(2*dt))*C)*Vt_x;
Gt_y = F_y + MandC*Rt_y + ((1/dt)*M - (1/(2*dt))*C)*Vt_y;

% solving the system
Rtplusdt_x = At\Gt_x;
Rtplusdt_y = At\Gt_y;

Rtminusdt_x = Rt_x;
Rtminusdt_y = Rt_y;
Rt_x = Rtplusdt_x;
Rt_y = Rtplusdt_y;

 % beginstate
mapToSpace(Rx,Ry);
AXIS([XMIN XMAX YMIN YMAX]);
Movie(frameNr) = getframe;
frameNr = frameNr + 1;
end