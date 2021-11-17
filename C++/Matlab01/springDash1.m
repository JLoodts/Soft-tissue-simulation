close all;
clear all;

m = 4.8e-7;
k = 1190;
c = 0.07; % 0.07 is nice
dt = 3e-7;
restLength = 0.003;

N = 1000;

t = zeros(N,1);
I = ones(size(t));

vel = 2;       % m/s double of the vel in SoftTissue
dir = 1;        % start with stretching

pos = restLength*I;
xs = 0*I;
vs = 0*I;
xd = 0*I;
vd = 0*I;

xs(1) = 0;
xd(1) = 0;
pos(1) = 0;
xs(2) = 0;
xd(2) = 0;
pos(2) = 0;
vs(1) = vel;

f = 0*I;
strain = 0*I;

strainMax = 5;

for i = 2:N-1
    
    t(i+1) = (i-1)*dt;
    pos(i+1) = pos(i) + vel*dir*dt;
    xs(i+1) = vel*2*dt-(k/c)*xs(i)*2*dt+xs(i-1);
    xd(i+1) = (k/c)*pos(i)*2*dt-(k/c)*xd(i)*2*dt+xd(i-1);
    
    vs(i) = (xs(i+1)-xs(i-1))/(2*dt);
    vd(i) = (xd(i+1)-xd(i-1))/(2*dt);

end

fs = k*xs;
fd = c*vd;

figure('Position',[0,0,500,700]);
subplot(3,1,1);
plot(t,xs,'-b'); hold on
plot(t,xd,'-r'); hold on
plot(t,xd+xs,'-g');
legend('xs','xd','x');
xlabel('time (s)');

subplot(3,1,2);
plot(t,vs,'-b'); hold on
plot(t,vd,'-r');
xlabel('time (s)');
% ylabel('velocity (m/s)');
legend('vs','vd');

subplot(3,1,3);
plot(t,fs,'ob'); hold on;
plot(t,fd,'-r');
xlabel('time (s)');
ylabel('force (N)');
legend('fs','fd');


figure('Position',[520,0,500,700]);
plot(xs+xd,fs,'ob'); hold on;
plot(xs+xd,fd,'-r');
xlabel('xs+xd (m)');
ylabel('force (N)');
legend('fs','fd');
%output = exactspringdash(0);