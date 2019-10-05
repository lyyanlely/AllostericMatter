function mm = BistableMat(x,pars) %nb,l0,flag,kw,a,b
% this function compute the energy and gradient for an elastic network
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

dim = size(x,1);
nnb = size(nb,1);

ii  = zeros(1,2*4*nnb);  %1:4*nnb, parallel, 4*nnb+1:2*4*nnb, perpen
jj  = zeros(1,2*4*nnb);
vv  = zeros(1,2*4*nnb);
kk  = zeros(1,2*nnb);   %1:nnb, parallel, nnb+1:2*nnb, perpen

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
    % parallel
    ii(4*(i-1)+1) = i;
    ii(4*(i-1)+2) = i;
    ii(4*(i-1)+3) = i;
    ii(4*(i-1)+4) = i;
    jj(4*(i-1)+1) = 2*n1-1;
    jj(4*(i-1)+2) = 2*n1;
    jj(4*(i-1)+3) = 2*n2-1;
    jj(4*(i-1)+4) = 2*n2;
    vv(4*(i-1)+1) = -dx(i)/dr(i);  % x
    vv(4*(i-1)+2) = -dy(i)/dr(i);
    vv(4*(i-1)+3) = dx(i)/dr(i);  % y
    vv(4*(i-1)+4) = dy(i)/dr(i);
    % perpendicular
    ii(4*nnb+4*(i-1)+1) = nnb+i;
    ii(4*nnb+4*(i-1)+2) = nnb+i;
    ii(4*nnb+4*(i-1)+3) = nnb+i;
    ii(4*nnb+4*(i-1)+4) = nnb+i;
    jj(4*nnb+4*(i-1)+1) = 2*n1-1;
    jj(4*nnb+4*(i-1)+2) = 2*n1;
    jj(4*nnb+4*(i-1)+3) = 2*n2-1;
    jj(4*nnb+4*(i-1)+4) = 2*n2;
    vv(4*nnb+4*(i-1)+1) = -dy(i)/dr(i);  % x
    vv(4*nnb+4*(i-1)+2) = dx(i)/dr(i);
    vv(4*nnb+4*(i-1)+3) = dy(i)/dr(i);  % y
    vv(4*nnb+4*(i-1)+4) = -dx(i)/dr(i);
    if flag(i) == 1
        e = dl(i)/a;
        kk(i) = kw*(3*e^2-1);  % effective stiffness
        kk(nnb+i) = kw*((e^2-1)*e+0.5*b)*a/dr(i);  % tension on the bond
    elseif flag(i) == -1
        e = dl(i)/a;
        kk(i) = kw*(3*e^2-1);  % effective stiffness
        kk(nnb+i) = kw*((e^2-1)*e-0.5*b)*a/dr(i);  % tension on the bond
    elseif flag(i) == 2
        e = dl(i)/a;
        kk(i) = kw*(3*e^2-1);  % effective stiffness
        kk(nnb+i) = kw*((e^2-1)*e)*a/dr(i);  % tension on the bond
    else
        kk(i) = 1;
        kk(nnb+i) = dl(i)/dr(i);
    end
end

ss = sparse(ii,jj,vv,2*nnb,dim);  % structure matrix
mm = ss'*diag(kk)*ss;