close all;
clear all;
% x-values
X = linspace(0,1,100);
% homogeneous Dirichlet boundary conditions
BC = [0,1,0;
      0,1,0];
% function-names e.g. F = @sin calls sin.m
F = @F;
G = @G;
R = @F;
Y = sturm(X,BC,F,G,R);
plot(X,Y);