function [f,g] = BisEnRot(x,pars) 
% this function compute the energy f and gradient force g for an elastic network
% pars.pos records positions of the fixed particles
% pars.nb records the connection of the free particles
% pars.flag records the weak connections
% pars.l0 records the rest length of springs between free particles
% pars.kw relative strength of weak springs
% pars.a records the rest length of bistable weak springs
% pars.b records the tilt of bistable weak springs
% pars.xb periodic boundary in x
% pars.yb periodic boundary in y
nb = pars.nb;

dx = x(2*nb(:,2)-1)-x(2*nb(:,1)-1);
dy = x(2*nb(:,2))-x(2*nb(:,1));

if isfield(pars,'xb')
    xbound = pars.xb;
    flg1 = dx>xbound/2;
    flg2 = dx<-xbound/2;
    dx(flg1) = dx(flg1)-xbound;
    dx(flg2) = dx(flg2)+xbound;
end
if isfield(pars,'yb')
    ybound = pars.yb;
    flg1 = dy>ybound/2;
    flg2 = dy<-ybound/2;
    dy(flg1) = dy(flg1)-ybound;
    dy(flg2) = dy(flg2)+ybound;
end

dr = sqrt(dx.^2+dy.^2);

dl = dr-pars.l0;

f = 0;
g = zeros(size(x));

nnb = size(nb,1);
if isfield(pars,'flag')
    flag= pars.flag;
    kw = pars.kw;
    a = pars.a;
    b = pars.b;
else
    flag= zeros(nnb,1);
end

for i = 1:nnb
    n1 = nb(i,1);
    n2 = nb(i,2);
    if flag(i) == 1
        e = dl(i)/a;
        f = f + 0.5*kw*a^2*(e^2*(0.5*e^2-1)+b*e);
        force = kw*((e^2-1)*dl(i)+0.5*a*b);
        g(2*n1-1) = g(2*n1-1)+force*(-dx(i))/dr(i);
        g(2*n1)   = g(2*n1)  +force*(-dy(i))/dr(i);
        g(2*n2-1) = g(2*n2-1)+force*dx(i)/dr(i);
        g(2*n2)   = g(2*n2)  +force*dy(i)/dr(i);
    elseif flag(i) == -1
        e = dl(i)/a;
        f = f + 0.5*kw*a^2*(e^2*(0.5*e^2-1)-b*e);
        force = kw*((e^2-1)*dl(i)-0.5*a*b);
        g(2*n1-1) = g(2*n1-1)+force*(-dx(i))/dr(i);
        g(2*n1)   = g(2*n1)  +force*(-dy(i))/dr(i);
        g(2*n2-1) = g(2*n2-1)+force*dx(i)/dr(i);
        g(2*n2)   = g(2*n2)  +force*dy(i)/dr(i);
    elseif flag(i) == 2
        e = dl(i)/a;
        f = f + 0.5*kw*a^2*e^2*(0.5*e^2-1);
        force = kw*(e^2-1)*dl(i);
        g(2*n1-1) = g(2*n1-1)+force*(-dx(i))/dr(i);
        g(2*n1)   = g(2*n1)  +force*(-dy(i))/dr(i);
        g(2*n2-1) = g(2*n2-1)+force*dx(i)/dr(i);
        g(2*n2)   = g(2*n2)  +force*dy(i)/dr(i);
    else
        f = f + 0.5*dl(i)^2;
        g(2*n1-1) = g(2*n1-1)+(-dx(i))/dr(i)*dl(i);
        g(2*n1)   = g(2*n1)  +(-dy(i))/dr(i)*dl(i);
        g(2*n2-1) = g(2*n2-1)+dx(i)/dr(i)*dl(i);
        g(2*n2)   = g(2*n2)  +dy(i)/dr(i)*dl(i);
    end
end

if isfield(pars,'mid')  % if there is a direction we don't want to relax along
    v = zeros(size(g));
    x1 = pars.mid-0.5-0.5/pars.L;
    x2 = pars.mid+0.5-0.5/pars.L;
    y0 = (1-1/pars.L)*sqrt(3)/4;
    id1 = find(x(1:2:end)<=pars.mid);
    id2 = find(x(1:2:end)>pars.mid);
    v(2*id1)   = x(2*id1-1)-x1;
    v(2*id1-1) = -(x(2*id1)-y0);
    v(2*id2)   = x(2*id2-1)-x2;
    v(2*id2-1) = -(x(2*id2)-y0);
    v = v/norm(v);
    g = g - v'*g*v;  % substract the force along the direction
end