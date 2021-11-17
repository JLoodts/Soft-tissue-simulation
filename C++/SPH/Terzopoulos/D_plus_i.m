function [R] = D_plus_i(i,R,h)

if (i==1)
    for m = 1:size(R,1)-1
        R(m,:) = (R(m+1,:) - R(m,:))./h;
    end
elseif (i==2)
    for n = 1:size(R,2)-1
        R(:,n) = (R(:,n+1) - R(:,n))./h;
    end
end
  