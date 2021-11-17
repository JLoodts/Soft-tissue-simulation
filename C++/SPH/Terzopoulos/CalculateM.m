function [Mass] = CalculateM(M,N)

length = M*N; 
Mass = diag( ones(1,length) ); 

