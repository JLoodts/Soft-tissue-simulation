function [V] = MatrixToVector(M)

% [V] = MatrixToVector(M) Transforms a matrix into a vector

[nrRows,nrColumns] = size(M);

V = [];
for i = 1:nrRows
    V = [V; M(i,:)'];
end