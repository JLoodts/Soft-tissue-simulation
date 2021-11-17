close all;
if (length(tout)~=length(x1sim))
    tout = tout(3:length(x1sim)+2);
end;
plot(tout,x1sim,'g');
hold on;
plot(tout,x2sim,'g');
hold on;

load forceDeformation.txt
t = forceDeformation(:,1);
x1 = forceDeformation(:,2);
v1 = forceDeformation(:,3);
a1 = forceDeformation(:,4);
x2 = forceDeformation(:,5);
v2 = forceDeformation(:,6);
a2 = forceDeformation(:,7);
clear forceDeformation;

plot(t,x1,'r');
hold on;
plot(t,x2,'r');
legend('Simulink x1','Simulink x2','simple scheme x1','simple scheme x2');
ylabel('displacement (m)');
xlabel('time (s)');
title('simulating a ystem with two masses');
% figure;
% subplot(2,1,1);
% plot(t,x1,'b-');
% hold on;
% plot(t,v1,'g-');
% hold on;
% plot(t,a1,'r-');
% subplot(2,1,2);
% plot(t,x2,'b-');
% hold on;
% plot(t,v2,'g-');
% hold on;
% plot(t,a2,'r-');
% xlabel('time (s)');
% title('simulating a double spring damper system');