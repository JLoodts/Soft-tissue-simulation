close all;
clear all;

load forceDeformation.txt
dt = forceDeformation(:,1);
err = forceDeformation(:,2);

clear forceDeformation;

figure;
err = err(2:length(err));
dt = dt(2:length(dt));
 semilogx(dt,err./0.001,'b-');
%loglog(dt,err./0.01,'b-');
xlabel('dt (s)');
ylabel('error (%)');
title('error in simulating a spring damper system');