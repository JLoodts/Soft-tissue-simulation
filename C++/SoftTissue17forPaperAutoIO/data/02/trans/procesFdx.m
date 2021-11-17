close all;
clear all;

load forceDeformation.txt
t = forceDeformation(:,1);
dx = forceDeformation(:,2);
f = forceDeformation(:,3);
vol = forceDeformation(:,4);
width = forceDeformation(:,5);
pressure = forceDeformation(:,6);
pressureMin = forceDeformation(:,7);
pressureMax = forceDeformation(:,8);
volCell = [];
nrCells = 42;
for i = 1:nrCells
  volCell = [volCell, forceDeformation(:,8+i)];
end
clear forceDeformation;

figure('Position',[0,0,500,700]);
subplot(3,1,1);
plot(dx,f,'-*');
title('simulated onion tissue');
%title('force-deformation curve for simulated onion tissue');
xlabel('strain (%)');
ylabel('stress (N/m²)');

subplot(3,1,2);
plot(t,dx,'+-');
%title('deformation with time for simulated onion tissue');
xlabel('time (s)');
ylabel('strain (%)');

subplot(3,1,3);
plot(t,vol,'-+'); hold on;
%title('deformation with time for simulated onion tissue');
xlabel('time (s)');
ylabel('volume (%)');%(m³)');

figure('Position',[520,0,500,700]);
subplot(3,1,1);
plot(dx,pressure,'b-*'); hold on;
plot(dx,pressureMin,'r-'); hold on;
plot(dx,pressureMax,'g-'); 
title('simulated onion tissue');
%title('force-deformation curve for simulated onion tissue');
xlabel('strain (%)');
ylabel('pressure (N/m²)');

subplot(3,1,2);
% [left, bottom, width, height]
for i = 1:nrCells
plot(t,volCell(:,i),'-'); hold on;
end
%title('deformation with time for simulated onion tissue');
xlabel('time (s)');
ylabel('volume (%)');%(m³)');


% figure('Position',[520,300,500,400]);
subplot(3,1,3);
% [left, bottom, width, height]
plot(dx,width,'-');
%title('deformation with time for simulated onion tissue');
xlabel('strain (%)');
ylabel('delta width (%)');%(m³)');

[P,S] = POLYFIT(dx,f,1);
P/1e6