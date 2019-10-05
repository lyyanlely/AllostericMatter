function q=nonlineartest(ff)

fname = '../Allostery/May23/Config_12_12_0.882_1.00_1.00_0.0500_0.0001_1_%03d.dat';
nr = 20;
nc = 5;
options.maxit = 10000;
options.prtlevel = 0;

xl    = 12;
yl    = 12;
dd    = 1;
dt    = 1.;
kweak = 0.0001;

np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
eps   = 0.2;   % distort size
warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl;
dim   = 2;
xbound = xl;
%% define perturb and target displacements
id  = (1:np);     %(xl-np)/2+ indices of the displaced particles
idd = dim*id(1)-1:dim*id(end);  %
it  = (yl-1)*xl+(1:nt);  %(yl-nt)/2+ indices of the target particles
itt = dim*it(1)-1:dim*it(end);
nid = setxor(1:dim*num,idd);
nit = setxor(1:dim*num,itt);
nidt= setxor(1:dim*num,[idd,itt]);
pert = zeros(dim*num,1);  % perturb displacement
targ = zeros(dim*num,1);  % target displacement
fixs = zeros(1,(np+nt)*2);
cf   = 0;
%% construct the embedded network
bondnb = zeros(3*num,2);
nnb    = 0;
bondwk = zeros(6*num,2);  % to all next neareat neighbors
nwk    = 0;
posx = zeros(1,num);
posy = zeros(1,num);
for n=1:num
    ny = ceil(n/xl);
    nx = n - (ny-1)*xl;
    posx(n) = nx-1+mod(ny-1,2)/2;
    posy(n) = sqrt(3)/2*ny - sqrt(3)/4;
    %
    if mod(ny,2)
        if mod(nx+mod(floor(ny/2),2),2)
            posx(n) = posx(n);
            posy(n) = posy(n);
        else
            posx(n) = posx(n) + eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    else
        if mod(nx+mod(ny/2-1,2),2)
            posx(n) = posx(n);
            posy(n) = posy(n) - eps;
        else
            posx(n) = posx(n) - eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    end
    %}
    %if nx<xl
    nnb = nnb+1;
    bondnb(nnb,1) = n;
    bondnb(nnb,2) = n+1-xl*(nx==xl);  % periodic in horizontal direction
    %end
    % fix the connections between the perturbing nodes and target nodes
    if ismember(n,id)&&ismember(n+1-xl*(nx==xl),id)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ismember(n,it)&&ismember(n+1-xl*(nx==xl),it)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ny<yl  % the top line is not periodic
        %if nx<xl
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl+(1-xl*(nx==xl))*mod(ny-1,2); %-num*(ny==linnum)
        %end
        %if nx>1
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl-(1-xl*(nx==1))*mod(ny,2);  %-num*(ny==linnum)
        %end
    end
    % weak bonds connect to the next nearest and next next nearest
    % next nearest neighbors
    if ny<yl
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl+1-xl*(nx==xl)+(1-xl*(nx==xl-1))*mod(ny-1,2);  %-num*(ny==linnum)
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl-1+xl*(nx==1)-(1-xl*(nx==2))*mod(ny,2);  %-num*(ny==linnum)
    end
    if ny<yl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl;
        % next-next nearest neighbors
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl+1-xl*(nx==xl);
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl-1+xl*(nx==1);
    end
    nwk = nwk+1;
    bondwk(nwk,1) = n;
    bondwk(nwk,2) = n+2-xl*(nx>=xl-1);
end
bondnb = bondnb(1:nnb,:);
bondwk = bondwk(1:nwk,:);
fixs   = fixs(1:cf);
pert(2*id(1)) = -1*dd;
targ(2*it(1)) = 1*dt;
for i = 2:np
    dx  = posx(id(i))-posx(id(i-1));
    dy  = posy(id(i))-posy(id(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    pert(2*id(i)-1) = pert(2*id(i-1)-1)*(dx^2-dy^2)+2*pert(2*id(i-1))*dx*dy;
    pert(2*id(i))   = pert(2*id(i-1))*(dy^2-dx^2)+2*pert(2*id(i-1)-1)*dx*dy;
end
for i = 2:nt
    dx  = posx(it(i))-posx(it(i-1));
    dy  = posy(it(i))-posy(it(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    targ(2*it(i)-1) = targ(2*it(i-1)-1)*(dx^2-dy^2)+2*targ(2*it(i-1))*dx*dy;
    targ(2*it(i))   = targ(2*it(i-1))*(dy^2-dx^2)+2*targ(2*it(i-1)-1)*dx*dy;
end
pos   = [posy;posx]';
pert1 = zeros(size(pert));
pert1(1:2:end) = pert(2:2:end);
pert1(2:2:end) = pert(1:2:end);
rid  = setdiff((1:num)',id);

%% Compute Mmatrix
bpx  = zeros(1,nnb);
bpy  = zeros(1,nnb);
delx = zeros(1,nnb);
dely = zeros(1,nnb);
delr = zeros(1,nnb);
idx = zeros(1,4*nnb);
jdx = zeros(1,4*nnb);
val = zeros(1,4*nnb);
for n = 1:nnb
    ii = bondnb(n,1);
    jj = bondnb(n,2);
    delx(n) = posx(ii) - posx(jj);
    dely(n) = posy(ii) - posy(jj);
    if delx(n)>xbound/2   % periodic
        delx(n) = delx(n) - xbound;
    elseif delx(n)<-xbound/2
        delx(n) = delx(n) + xbound;
    end
    bpx(n) = posx(jj)+delx(n)/2;
    bpy(n) = posy(jj)+dely(n)/2;
    delr(n) = sqrt(delx(n)^2+dely(n)^2);
    idx(4*(n-1)+1) = n;
    jdx(4*(n-1)+1) = 2*(ii-1)+1;
    val(4*(n-1)+1) = delx(n)/delr(n);
    idx(4*(n-1)+2) = n;
    jdx(4*(n-1)+2) = 2*(ii-1)+2;
    val(4*(n-1)+2) = dely(n)/delr(n);
    idx(4*(n-1)+3) = n;
    jdx(4*(n-1)+3) = 2*(jj-1)+1;
    val(4*(n-1)+3) = -delx(n)/delr(n);
    idx(4*(n-1)+4) = n;
    jdx(4*(n-1)+4) = 2*(jj-1)+2;
    val(4*(n-1)+4) = -dely(n)/delr(n);
end
Smatrix = sparse(idx,jdx,val,nnb,dim*num);

% of the weak network
bpxw = zeros(1,nnb);
bpyw = zeros(1,nnb);
delxw = zeros(1,nwk);
delyw = zeros(1,nwk);
delrw = zeros(1,nwk);
idxw = zeros(1,4*nwk);
jdxw = zeros(1,4*nwk);
valw = zeros(1,4*nwk);
for n = 1:nwk
    ii = bondwk(n,1);
    jj = bondwk(n,2);
    delxw(n) = posx(ii) - posx(jj);
    delyw(n) = posy(ii) - posy(jj);
    if delxw(n)>xbound/2
        delxw(n) = delxw(n) - xbound;
    elseif delxw(n)<-xbound/2
        delxw(n) = delxw(n) + xbound;
    end
    bpxw(n) = posx(ii)+delxw(n)/2;
    bpyw(n) = posy(ii)+delyw(n)/2;
    delrw(n) = sqrt(delxw(n)^2+delyw(n)^2);
    idxw(4*(n-1)+1) = n;
    jdxw(4*(n-1)+1) = 2*(ii-1)+1;
    valw(4*(n-1)+1) = delxw(n)/delrw(n);
    idxw(4*(n-1)+2) = n;
    jdxw(4*(n-1)+2) = 2*(ii-1)+2;
    valw(4*(n-1)+2) = delyw(n)/delrw(n);
    idxw(4*(n-1)+3) = n;
    jdxw(4*(n-1)+3) = 2*(jj-1)+1;
    valw(4*(n-1)+3) = -delxw(n)/delrw(n);
    idxw(4*(n-1)+4) = n;
    jdxw(4*(n-1)+4) = 2*(jj-1)+2;
    valw(4*(n-1)+4) = -delyw(n)/delrw(n);
end
SMatW = sparse(idxw,jdxw,valw,nwk,dim*num);
MmatW = SMatW'*SMatW;

q = zeros(length(ff),nr*nc);
for r = 1:20
    filename = sprintf(fname,r);
    config   = dlmread(filename);
    ll = size(config,1);
    lc = ll/2/nc;
    for ii = 1:nc
        stat = config(ll/2+lc*ii,:);
%% linear response
pMmat  = Smatrix'*diag(stat)*Smatrix+kweak*MmatW;
pQmat   = zeros(dim*num);
for i=idd
    pQmat(i,i) = 1;
end
pQmat(:,nid) = -pMmat(:,nid);
fext = pMmat*(pert);
disp = linsolve(pQmat,fext);  %pQmat\fext; %
if any(isnan(disp)|abs(disp)>1e5)
    disp = pinv(pQmat)*fext; %
end
disp(idd) = pert(idd);
xid = false(dim*num,1);
yid = false(dim*num,1);
xid(1:2:dim*num) = true;
yid(2:2:dim*num) = true; %mod(nid,2)==0;
dx = mean(disp((xid))); %nid
dy = mean(disp((yid)));
displ = disp-dx*xid-dy*yid; %nid

%% conjugate gradient
pars = initialpars(num,pos,bondnb(stat>0,:),id,xbound);
options.x0 = pars.x0;
pars.nvar  = dim*(num-length(id));
pars.fgname= 'networkcg';
for j = 1:length(ff)
options.normtol=1e-6*ff(j);
pars.pos   = pars.pos+pert1(idd)*ff(j);
[x,~,~,~,~]=nlcg(pars,options);
xx = reshape(x,2,[]);
px = xx(2,:);
py = xx(1,:);
dpx = px-posx(rid);
dpy = py-posy(rid);
dpnl  = zeros(size(displ));
dpnl(idd) = pert(idd)*ff(j);
dpnl(nid(1:2:end)) = dpx';
dpnl(nid(2:2:end)) = dpy';
dx = mean([dpx,ff(j)*pert(idd(1:2:end))']);
dy = mean([dpy,ff(j)*pert(idd(2:2:end))']);
dpnl  = dpnl-dx*xid-dy*yid;

q(j,nc*(r-1)+ii) = dpnl(itt)'*displ(itt)/norm(dpnl(itt))/norm(displ(itt));
end
    end
end