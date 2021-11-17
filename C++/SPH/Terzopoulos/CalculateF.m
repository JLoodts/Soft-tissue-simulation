function [Fx,Fy] = CalculateF(M,N)

length = M*N; 
Fx = -1*ones(length,1); 
Fy = -9.81*ones(length,1); 
% Fx(1,1) = 2;
% Fy(1,1) = Fy(1,1)-2;