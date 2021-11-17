clear all;
close all;
load forceDeformation.txt
t = forceDeformation(:,1);
x1 = forceDeformation(:,2);
v1 = forceDeformation(:,3);
a1 = forceDeformation(:,4);
x2 = forceDeformation(:,5);
v2 = forceDeformation(:,6);
a2 = forceDeformation(:,7);
f1 = forceDeformation(:,8);
f2 = forceDeformation(:,9);

clear forceDeformation;

% plot(t,x1,'r');
% hold on;
% plot(t,x2,'b');
% ylabel('displacement (m)');
% xlabel('time (s)');
% title('simulating a system with two masses');
plot(x2,f1,'r');
hold on;
plot(x2,f2,'b');
xlabel('displacement (m)');
ylabel('force (N)');
title('simulating a system with two masses');

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