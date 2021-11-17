close all;
clear all;
syms h1x h1y h2x h2y h3x h3y
a = sqrt((h2x-h1x)^2 + (h2y-h1y)^2);
b = sqrt((h3x-h2x)^2 + (h3y-h2y)^2);
c = sqrt((h1x-h3x)^2 + (h1y-h3y)^2);
s = ((a+b+c)/2);
surface = sqrt(s*(s-a)*(s-b)*(s-c));
expand(surface)