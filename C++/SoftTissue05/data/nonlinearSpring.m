clear all;
close all;
t = linspace(0,0.001,1000);
% % x = 10000*t.^1./(ones(size(t))+10000*t.^1);
% % parabool = (2.5e-5-(100*(t-0.0005*ones(size(t))).^2));
goniometric = t.*(cos(1000*t));
%parabool = -0.000025+1000*(t-0.0005*ones(size(t))).^2;
x = t + goniometric;

plot(t,x);
hold on;
plot(t,goniometric,'r');
title('goniometric spring');
legend('force','extra term');
xlabel('displacement');