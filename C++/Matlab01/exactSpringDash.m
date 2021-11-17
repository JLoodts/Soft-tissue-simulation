
function [result] = exactspringdash(input)

close all;
clear all;

m = 4.8e-7;
k = 1190;
c = 0.07; % 0.01 is nice
% dt = 3e-7;
restLength = 0.003;

N = 1000;

% t = zeros(N,1);
t = linspace(0,3e-4,N); % s
dt = t(2)-t(1)
I = ones(size(t));

vel = 2;       % m/s double of the vel in SoftTissue

xs = 0*I;
vs = 0*I;
xd = 0*I;
vd = 0*I;

xs0 = 0;
xd0 = 0;
fs = 0*I;
fd = 0*I;

xd = -(vel*c/k)*I + vel*t + (xd0+(vel*c/k))*exp(-(k/c)*t);
vd = vel*I - k/c*(xd0+(vel*c/k))*exp(-(k/c)*t);
xs = (vel*c/k)*I + (xs0-(vel*c/k))*exp(-(k/c)*t);
vs = -k/c*(xs0-(vel*c/k))*exp(-(k/c)*t);
fs = k*xs;
fd = c*vd;

figure('Position',[0,0,500,700]);
subplot(3,1,1);
plot(t,xs,'-b'); hold on;
plot(t,xd,'-r'); hold on;
plot(t,xs+xd,'-g');
legend('xs','xd','x');
xlabel('time (s)');

subplot(3,1,2);
plot(t,vs,'-b'); hold on
plot(t,vd,'-r');
xlabel('time (s)');
% ylabel('velocity (m/s)');
legend('vs','vd');

subplot(3,1,3);
plot(t,fs,'+b'); hold on;
plot(t,fd,'-r');
xlabel('time (s)');
ylabel('force (N)');
legend('fs','fd');


figure('Position',[520,0,500,700]);
plot(xs+xd,fs,'+b'); hold on;
plot(xs+xd,fd,'-r');
xlabel('xs+xd (m)');
ylabel('force (N)');
legend('fs','fd');
result = 0;