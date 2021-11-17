close all;
clear all;

load forceDeformation.txt
t = forceDeformation(:,1);
dx = forceDeformation(:,2);
f = forceDeformation(:,3);
clear forceDeformation;

subplot(2,1,1);
plot(dx,f,'-');
title('simulated onion tissue');
%title('force-deformation curve for simulated onion tissue');
xlabel('displacement (m)');
ylabel('force (N)');

subplot(2,1,2);
plot(t,dx,'-');
%title('deformation with time for simulated onion tissue');
xlabel('time (s)');
ylabel('displacement (m)');
