function [R] = D_plus_ij(i,j,R)

R_temp = R;
if ((i==1)&(j==1))
    h1_pow2 = (R(2,1)-R(1,1))^2;
    for m = 2:size(R,1)-1
        R(m,:) = (R_temp(m+1,:) - 2.*R_temp(m,:) + R_temp(m-1,:))./h1_pow2; 
    end
elseif ((i==2)&(j==2))
    h2_pow2 = (R(1,2)-R(1,1))^2;
    for n = 2:size(R,2)-1
        R(:,n) = (R_temp(:,n+1) - 2.*R_temp(:,n) + R_temp(:,n-1))./h2_pow2; 
    end
elseif (((i==1)&(j==2))|((i==2)&(j==1)))
    h1_h2 = (R(2,1)-R(1,1))*(R(1,2)-R(1,1));
    for n = 2:size(R,2)-1
        for m = 2:size(R,1)-1
            R(m,n) = (R_temp(m+1,n+1) - R_temp(m+1,n) - R_temp(m,n+1) + R_temp(m,n))./h1_h2; 
       end
    end
end