function [surface] = calculateTriangleSurface2(xpos, ypos)
h1x = xpos(1); h1y = ypos(1);
h2x = xpos(2); h2y = ypos(2);
h3x = xpos(3); h3y = ypos(3);
surface = (1/4*h2x^2*h1y^2+1/4*h2x^2*h3y^2-1/2*h2x*h1x*h2y*h1y+1/2*h1y*h3y*h3x*h2x+1/2*h1x*h3x*h2y*h1y+1/2*h2x*h1x*h3y*h2y+1/2*h2x*h1x*h1y*h3y+1/2*h2y*h1y*h3x*h2x-1/2*h3x*h2x*h3y*h2y+1/2*h3y*h2y*h1x*h3x-1/2*h1x*h3x*h1y*h3y+1/4*h1x^2*h3y^2-1/2*h1y*h3y*h2x^2-1/2*h3x*h2x*h1y^2-1/2*h1x*h3x*h2y^2-1/2*h2x*h1x*h3y^2-1/2*h2y*h1y*h3x^2-1/2*h3y*h2y*h1x^2+1/4*h3x^2*h2y^2+1/4*h3x^2*h1y^2+1/4*h2y^2*h1x^2)^(1/2);