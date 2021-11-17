function [C] = CalculateC(M,N)

length = M*N; 
C = diag( 0.2*ones(1,length) ); 
