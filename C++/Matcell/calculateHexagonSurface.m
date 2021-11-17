function [surface] = calculateHexagonSurface(pos)

% pos = [x1,y1,x2,y2...]
surface = 0;
for i = 2:5
    surface = surface + calculateTriangleSurface(...
        [pos(1),    pos(2*i-1), pos(2*(i+1)-1)  ],...
        [pos(2),    pos(2*i),   pos(2*(i+1))    ]);
end