function [f,g] = BistableE(x,pars) 
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
% pars.v direction to be removed from relaxation
% pars.bid indecies of binding sites
% pars.k0 strength of constraints on binding sites
% pars.th0 reference angle of constraint
% pars.r0  reference distance of binding sites
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
    if flag(i) == 0
        f = f + 0.5*dl(i)^2;
        g(2*n1-1) = g(2*n1-1)+(-dx(i))/dr(i)*dl(i);
        g(2*n1)   = g(2*n1)  +(-dy(i))/dr(i)*dl(i);
        g(2*n2-1) = g(2*n2-1)+dx(i)/dr(i)*dl(i);
        g(2*n2)   = g(2*n2)  +dy(i)/dr(i)*dl(i);
    else
        if isfield(pars,'p')
            e = dr(i)/a;
            f = f + pars.p(i,1)*e^4+pars.p(i,2)*e^3+pars.p(i,3)*e^2+pars.p(i,4)*e+pars.p(i,5);
            force = (4*pars.p(i,1)*e^3+3*pars.p(i,2)*e^2+2*pars.p(i,3)*e+pars.p(i,4))/a;
            g(2*n1-1) = g(2*n1-1)+force*(-dx(i))/dr(i);
            g(2*n1)   = g(2*n1)  +force*(-dy(i))/dr(i);
            g(2*n2-1) = g(2*n2-1)+force*dx(i)/dr(i);
            g(2*n2)   = g(2*n2)  +force*dy(i)/dr(i);
        else
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
            end
        end
    end
end


if isfield(pars,'bid')
    n1 = pars.bid(1);
    n2 = pars.bid(2);
    dx = x(2*n2-1)-x(2*n1-1);
    dy = x(2*n2)-x(2*n1);
    dr = sqrt(dx^2+dy^2);
    th = atan(dy/dx);
    if dx<0
        th = th+pi;
    end
    fth = pars.k0*(th-pars.th0)*pars.l0^2;
    fr  = pars.k0*(dr-pars.r0);
    f = f + 0.5*pars.k0*(th-pars.th0)^2*pars.l0^2;
    f = f + 0.5*pars.k0*(dr-pars.r0)^2;
    g(2*n1-1) = g(2*n1-1)+fr*(-dx)/dr;
    g(2*n1)   = g(2*n1)  +fr*(-dy)/dr;
    g(2*n2-1) = g(2*n2-1)+fr*dx/dr;
    g(2*n2)   = g(2*n2)  +fr*dy/dr;
    g(2*n1-1) = g(2*n1-1)+fth*dy/dr^2;
    g(2*n1)   = g(2*n1)  +fth*(-dx)/dr^2;
    g(2*n2-1) = g(2*n2-1)+fth*(-dy)/dr^2;
    g(2*n2)   = g(2*n2)  +fth*dx/dr^2;
end

if isfield(pars,'bit')
    n1 = pars.bit(1);
    n2 = pars.bit(2);
    dx = x(2*n2-1)-x(2*n1-1);
    dy = x(2*n2)-x(2*n1);
    dr = sqrt(dx^2+dy^2);
    th = atan(dy/dx);
    if dx<0
        th = th+pi;
    end
    fth = pars.k0*(th-pars.th1)*pars.l0^2;
    fr  = pars.k0*(dr-pars.r1);
    f = f + 0.5*pars.k0*(th-pars.th1)^2*pars.l0^2;
    f = f + 0.5*pars.k0*(dr-pars.r1)^2;
    g(2*n1-1) = g(2*n1-1)+fr*(-dx)/dr;
    g(2*n1)   = g(2*n1)  +fr*(-dy)/dr;
    g(2*n2-1) = g(2*n2-1)+fr*dx/dr;
    g(2*n2)   = g(2*n2)  +fr*dy/dr;
    g(2*n1-1) = g(2*n1-1)+fth*dy/dr^2;
    g(2*n1)   = g(2*n1)  +fth*(-dx)/dr^2;
    g(2*n2-1) = g(2*n2-1)+fth*(-dy)/dr^2;
    g(2*n2)   = g(2*n2)  +fth*dx/dr^2;
end

if isfield(pars,'v')  % if there is a direction we don't want to relax along
    g = g - pars.v'*g*pars.v;  % substract the force along the direction
end
%{
if isfield(pars,'bp0')
    dg = -pars.k0*(pars.bp0-x(pars.bid));
    dg(1:2:end) = 0;
    dg(2:2:end) = dg(2:2:end)-mean(dg(2:2:end));
    g(pars.bid) = g(pars.bid)+dg;
    %f = f+0.5*pars.k0*sum((pars.bp0-x(pars.bid)).^2);
end
%}