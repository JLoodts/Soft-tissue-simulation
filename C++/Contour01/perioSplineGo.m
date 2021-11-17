clear all;
close all;

load clickResult.mat;
figure;
plot(x,y,'b+');
hold on;

f = [x;y];
n = length(x);
clear x y;
x = linspace(0,1,n+1);
t = linspace(0,1,100);
y = periospline(x,f,t);

plot(y(1,:),y(2,:),'r-');
axis equal;
figure;
plot(t(1:size(y,2)),y,'b-');