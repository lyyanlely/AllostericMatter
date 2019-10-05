function [x,v] = shearsquare(x,theta,xl,x0)
% shear (rotate) a block cut at x0 (center at 0.5+x0) by angle theta

if nargin<4
    x0 = 0.5;
end
xtemp = x;
v = zeros(size(x));

x1 = x0-0.5-0.5/xl;
x2 = x0+0.5-0.5/xl;
y0 = (1-1/xl)*sqrt(3)/4;

id1 = find(x(1:2:end)<=x0);
id2 = find(x(1:2:end)>x0);

x(2*id1-1) = x1+cos(theta)*(xtemp(2*id1-1)-x1)-sin(theta)*(xtemp(2*id1)-y0);
v(2*id1)   = cos(theta)*(xtemp(2*id1-1)-x1)-sin(theta)*(xtemp(2*id1)-y0);
x(2*id1)   = y0+sin(theta)*(xtemp(2*id1-1)-x1)+cos(theta)*(xtemp(2*id1)-y0);
v(2*id1-1) = -(sin(theta)*(xtemp(2*id1-1)-x1)+cos(theta)*(xtemp(2*id1)-y0));

x(2*id2-1) = x2+cos(theta)*(xtemp(2*id2-1)-x2)-sin(theta)*(xtemp(2*id2)-y0);
v(2*id2)   = cos(theta)*(xtemp(2*id2-1)-x2)-sin(theta)*(xtemp(2*id2)-y0);
x(2*id2)   = y0+sin(theta)*(xtemp(2*id2-1)-x2)+cos(theta)*(xtemp(2*id2)-y0);
v(2*id2-1) = -(sin(theta)*(xtemp(2*id2-1)-x2)+cos(theta)*(xtemp(2*id2)-y0));