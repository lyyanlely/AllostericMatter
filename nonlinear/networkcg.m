function [f,g] = networkcg(x,pars)
% this function compute the energy and gradient for an elastic network
% pars.pos records positions of the fixed particles
% pars.nb records the connection of the free particles
% pars.ext records the connection between the fixed and free particles
% pars.nbl records the rest length of springs between free particles
% pars.exl records the rest length of springs between fixed and free
% particles

nnb = size(pars.nb,1);
nex = size(pars.ext,1);

L = pars.L;

f = 0;
g = zeros(size(x));

for i = 1:nnb
    n1   = pars.nb(i,1);
    n2   = pars.nb(i,2);
    delx = x(2*n2-1)-x(2*n1-1);
    dely = x(2*n2)  -x(2*n1);
    if dely>.5*L
        dely = dely-L;
    elseif dely<-.5*L
        dely = dely+L;
    end
    delr = sqrt(delx^2+dely^2);
    difr = delr-pars.nbl(i);
    f = f + 0.5*difr^2;
    g(2*n1-1) = g(2*n1-1)+(-delx)/delr*difr;
    g(2*n1)   = g(2*n1)  +(-dely)/delr*difr;
    g(2*n2-1) = g(2*n2-1)+delx/delr*difr;
    g(2*n2)   = g(2*n2)  +dely/delr*difr;
end

for i = 1:nex
    n0   = pars.ext(i,1);
    n1   = pars.ext(i,2);
    delx = x(2*n1-1)-pars.pos(2*n0-1);
    dely = x(2*n1)  -pars.pos(2*n0);
    if dely>.5*L
        dely = dely-L;
    elseif dely<-.5*L
        dely = dely+L;
    end
    delr = sqrt(delx^2+dely^2);
    difr = delr-pars.exl(i);
    f = f + 0.5*difr^2;
    g(2*n1-1) = g(2*n1-1)+delx/delr*difr;
    g(2*n1)   = g(2*n1)  +dely/delr*difr;
end