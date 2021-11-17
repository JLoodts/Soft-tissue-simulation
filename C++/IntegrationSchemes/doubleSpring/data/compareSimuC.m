close all;

load forceDeformation.txt
t = forceDeformation(:,1);
x = forceDeformation(:,2);
v = forceDeformation(:,3);
a = forceDeformation(:,4);
x_e = forceDeformation(:,5);
% v_e = forceDeformation(:,6);
% a_e = forceDeformation(:,7);
clear forceDeformation;

xerr = xsim - x_e;
plot(t,xerr,'r+-');
