function [K] = CalculateK(Rx,Ry,M,N,h1,h2)

% material parameters
nu_11 = 50*ones(M,N);
nu_12 = 10*ones(M,N);
nu_21 = nu_12;
nu_22 = 20*ones(M,N);


lengthK = M*N; 
K = spalloc( lengthK, lengthK, 7*lengthK ); % stiffness matrix
    % maximum memory requirement is k*lengthK because the stencil has k elements

% the degrees of freedom are stored in the column vector R, for which holds:
%    1)  K * Ry = Fy and K * Ry = Fy
%    2)  R(k) = r(i,j) with r(i,j) the mapping of point (i,j)



h1_pow2 = h1^2;
h1_pow4 = h1^4;
h2_pow2 = h2^2;
h2_pow4 = h2^4;
h1_h2   = h1*h2;
h1_h2_pow2 = h1_h2^2;

% given R0, calculate G_ij and B_ij
D_1_plusX = D_plus_i(1,Rx,h1); D_1_plusY = D_plus_i(1,Ry,h1);
D_2_plusX = D_plus_i(2,Rx,h2); D_2_plusY = D_plus_i(2,Ry,h2);
G_11_zero = D_1_plusX.*D_1_plusX + D_1_plusY.*D_1_plusY;
G_12_zero = D_1_plusX.*D_2_plusX + D_1_plusY.*D_2_plusY;
G_21_zero = G_12_zero;
G_22_zero = D_2_plusX.*D_2_plusX + D_2_plusY.*D_2_plusY;

alfa_11 = (D_1_plusX.*D_1_plusX + D_1_plusY.*D_1_plusY - G_11_zero).*nu_11;
alfa_12 = (D_1_plusX.*D_2_plusX + D_1_plusY.*D_2_plusY - G_12_zero).*nu_12;
alfa_21 = alfa_12;
alfa_22 = (D_2_plusX.*D_2_plusX + D_2_plusY.*D_2_plusY - G_22_zero).*nu_22;

% filling K
 for m = 2:M-1 
    for n = 2:N-1 

       k = N*(m-1) + n; % rownumber of the corresponding degree of freedom in K

       a = alfa_12(m-1,  n)/h1_h2;
       b = alfa_12(  m,  n)/h1_h2;
       c = alfa_22(  m,  n)/h2_pow2;
       d = alfa_11(m-1,  n)/h1_pow2;
       e = alfa_11(  m,  n)/h1_pow2;
       f = alfa_11(m-1,  n)/h1_pow2;
       g = alfa_22(  m,n-1)/h2_pow2;
       h = alfa_11(  m,  n)/h1_pow2;
       i = alfa_21(  m,n-1)/h1_h2;
           
       K(k,k-1-N) = a; % refers to Element(m-1,n+1)
       K(k,k-N)   = -b-c;
       K(k,k-1)   = -a-d;
       K(k,k)     = 2*b+c+e+f+g;
       K(k,k+1)   = -b-h;
       K(k,k+N)   = -g-i;
       K(k,k+1+N) = i;
       
    end 
 end 