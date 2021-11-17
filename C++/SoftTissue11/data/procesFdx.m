close all;
clear all;

load forceDeformation.txt
t = forceDeformation(:,1);
dx = forceDeformation(:,2);
f = forceDeformation(:,3);
vol = forceDeformation(:,4);
clear forceDeformation;

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
plot(t,vol,'+-');
%title('deformation with time for simulated onion tissue');
xlabel('time (s)');
ylabel('volume (%)');%(m³)');