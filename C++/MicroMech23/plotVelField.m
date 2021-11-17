% plotVelField


clear;
close all;
load outputX.txt;
load outputY.txt;
%z = sqrt(outputX.^2 + outputY.^2);
zx = outputX;
zy = outputY;
clear outputX outputY;
[dimY,dimX] = size(zx);
x = [1:dimX];%linspace(0,0.08,dimX);
y = [1:dimY];%linspace(0,0.0101,dimY);

%y = outputY;


%surfl(x,y,z); 
pcolor(x,y,zx);
shading INTERP;  % FACETED FLAT INTERP 
axis equal;
axis([min(x),max(x),min(y),max(y)]);
title('snelheidsprofiel in x-richting');
colorbar('horiz');
colormap(jet);
figure;
pcolor(x,y,zy);
shading INTERP;  % FACETED FLAT INTERP 
axis equal;
axis([min(x),max(x),min(y),max(y)]);
title('snelheidsprofiel in y-richting');
colorbar('horiz');
colormap(jet);