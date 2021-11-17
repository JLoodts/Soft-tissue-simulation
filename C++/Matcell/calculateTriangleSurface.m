function [surface] = calculateTriangleSurface(xpos, ypos)

a = sqrt((xpos(2)-xpos(1))^2 + (ypos(2)-ypos(1))^2);
b = sqrt((xpos(3)-xpos(2))^2 + (ypos(3)-ypos(2))^2);
c = sqrt((xpos(1)-xpos(3))^2 + (ypos(1)-ypos(3))^2);
s = 0.5*(a+b+c);
surface = sqrt(s*(s-a)*(s-b)*(s-c));