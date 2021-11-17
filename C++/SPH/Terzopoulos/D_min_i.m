function [R] = D_min_i(i,R)

if (i==1)
    h1 = (R(2,1)-R(1,1));
    for m = length(R,1):2
        R(m,:) = R(m,:) - R(m-1,:))./h1;
    end
elseif (i==2)
    h2 = (R(1,2)-R(1,1));
    for n = length(R,2):2
        R(:,n) = R(:,n) - R(:,n-1))./h2;
    end
end