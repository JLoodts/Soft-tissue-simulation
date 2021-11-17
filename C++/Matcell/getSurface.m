function [surface] = getSurface(Node)

n_max = 3;
A = [1, 2, 3];
B = [2, 3, 4];
C = [5, 5, 5];
surface = 0;
for n = 1:n_max
    surface = surface + calculateTriangleSurface(...
        [Node(A(n),2),    Node(B(n),2),   Node(C(n),2)  ],...
        [Node(A(n),3),    Node(B(n),3),   Node(C(n),3)  ]);
end