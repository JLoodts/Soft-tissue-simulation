nrCoo = length(x);
rightForm = [];
for i=1:nrCoo
    rightForm = [rightForm, x(i), y(i)];
end
    
save cornerPoints.txt rightForm -ASCII