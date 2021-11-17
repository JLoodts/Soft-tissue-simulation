function [M] = VectorToMatrix(V,nrRows,nrColumns)

% [M] = VectorToMatrix(V,nrRows,nrColumns) Transforms a vector into a matrix

M = zeros(nrRows,nrColumns);

for i = 1:nrRows
    beginPoint = 1+(i-1)*nrColumns; endPoint = beginPoint+nrColumns-1; 
    M(i,:) = V([beginPoint:endPoint])';
end