function allostery1(replica, xl, yl, spf, beta_ini)
% similar function as allostery, able to choose cooperative cost or geometric loss by changing parameter "allos"
%replica = 1;
allos = 0;  % 0 for cooperative 1 for geometric
%xl    = 12;
%yl    = 12;
%spf   = 0.83;
%nsp   = 385;%
%i0    = 3;
dd    = 1;
dt    = 1.;
kweak = 0.0001;
%beta_ini = 0.0001;
drbin = 0.2;
topnb = 36;

lmem  = 32;
%dmem  = 'uint32';

np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
temprange = 1;
iteration = 1000;
paraiter  = 1000;
interv = 100;
nrec   = round(iteration/interv);
prec   = round(paraiter/interv);
npara  = round(iteration/paraiter);
eps   = 0.2;   % distort size
%confint = 100;
%slow = 10;
warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl;
dim   = 2;
xbound = xl;
ybound = yl*sqrt(3)/2;
nmax  = 3*num-2*xl;
nflip = ones(1,temprange);
nf    = ones(1,temprange);
beta  = 1/beta_ini;
for tt = 1:temprange
    nflip(tt) = temprange+1-tt;
end

%% define perturb and target displacements
id  = (xl-np)/2+(1:np);     % indices of the displaced particles
idd = dim*id(1)-1:dim*id(end);  %
it  = (xl-nt)/2+(yl-1)*xl+(1:nt);  % indices of the target particles
itt = dim*it(1)-1:dim*it(end);
idt = [idd,itt];
nid = setxor(1:dim*num,idd);
nit = setxor(1:dim*num,itt);
nidt= setxor(1:dim*num,[idd,itt]);
ndd = 2*dim+1:dim*num;
pts  = zeros(dim*num,dim);  % translation on perturbing nodes
tts  = zeros(dim*num,dim);  % translation on target nodes
pert = zeros(dim*num,1);  % perturb displacement
targ = zeros(dim*num,1);  % target displacement
ptx  = zeros(dim*num,1);  % translation in x direction on perturbing nodes
pty  = zeros(dim*num,1);  % translation in y direction on perturbing
prot = zeros(dim*num,1);  % rotation on perturbing
ttx  = zeros(dim*num,1);  % translation in x direction on target nodes
tty  = zeros(dim*num,1);  % translation in y direction on target
trot = zeros(dim*num,1);  % rotation on target
tx   = zeros(dim*num,1);  % translation in x direction on target nodes
ty   = zeros(dim*num,1);  % translation in y direction on target
rot  = zeros(dim*num,1);  % rotation on target
fixs = zeros(1,(np+nt)*2);  % labels of springs to be fixed
cf   = 0;
%% construct the embedded network
bondnb = zeros(3*num,2);
nnb    = 0;
bondwk = zeros(6*num,2);  % to all next neareat neighbors
nwk    = 0;
posx = zeros(1,num);
posy = zeros(1,num);
for nn=1:num
    ny = ceil(nn/xl);
    nx = nn - (ny-1)*xl;
    n = nn;
    %nx = mod(nx+7,xl)+1;
    %n  = (ny-1)*xl+nx;
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
    %theta = rand(1);
    %posx(n) = posx(n) + eps*cos(theta);
    %posy(n) = posy(n) + eps*sin(theta);
    if nx<xl
    nnb = nnb+1;
    bondnb(nnb,1) = n;
    bondnb(nnb,2) = n+1-xl*(nx==xl);  % periodic in horizontal direction
    end
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
        if nx<xl ||mod(ny,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl+(1-xl*(nx==xl))*mod(ny-1,2); %-num*(ny==linnum)
        end
        if nx>1 ||mod(ny+1,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl-(1-xl*(nx==1))*mod(ny,2);  %-num*(ny==linnum)
        end
    end
    % weak bonds connect to the next nearest and next next nearest
    % next nearest neighbors
    if ny<yl
        if nx<xl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl+1-xl*(nx==xl)+(1-xl*(nx==xl-1))*mod(ny-1,2);  %-num*(ny==linnum)
        end
        if nx>2
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl-1+xl*(nx==1)-(1-xl*(nx==2))*mod(ny,2);  %-num*(ny==linnum)
        end
    end
    if ny<yl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl;
        % next-next nearest neighbors
        if nx<xl
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl+1-xl*(nx==xl);
        end
        if nx>1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl-1+xl*(nx==1);
        end
    end
    if nx<xl-1
    nwk = nwk+1;
    bondwk(nwk,1) = n;
    bondwk(nwk,2) = n+2-xl*(nx>=xl-1);
    end
end
bondnb = bondnb(1:nnb,:);
bondwk = bondwk(1:nwk,:);
fixs   = fixs(1:cf);
%% initialize the perturbing and target displacements
pert(2*id(1)) = 1*dd;
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
% translational and rotational parts at the stimulus site and target site
pcen = [mean(posx(id)),mean(posy(id))]; % center of perturbing nodes
tcen = [mean(posx(it)),mean(posy(it))]; % center of targeting nodes
ptc  = [mean(posx([id,it])),mean(posy([id,it]))]; % center of stimulus and target
ptx(idd(1:2:end)) = 1;
pty(idd(2:2:end)) = 1;
ttx(itt(1:2:end)) = 1;
tty(itt(2:2:end)) = 1;
tx(idd(1:2:end))  = 1;
tx(itt(1:2:end))  = 1;
ty(idd(2:2:end))  = 1;
ty(itt(2:2:end))  = 1;
for i = 1:np
    dx = posx(id(i)) - pcen(1);
    dy = posy(id(i)) - pcen(2);
    prot(2*id(i)-1) = -dy;
    prot(2*id(i))   = dx;
    dx = posx(id(i)) - ptc(1);
    dy = posy(id(i)) - ptc(2);
    rot(2*id(i)-1)  = -dy;
    rot(2*id(i))    = dx;
end
for i = 1:nt
    dx = posx(it(i)) - tcen(1);
    dy = posy(it(i)) - tcen(2);
    trot(2*it(i)-1) = -dy;
    trot(2*it(i))   = dx;
    dx = posx(it(i)) - ptc(1);
    dy = posy(it(i)) - ptc(2);
    rot(2*it(i)-1)  = -dy;
    rot(2*it(i))    = dx;
end
ptx  = ptx/norm(ptx);
pty  = pty/norm(pty);
prot = prot/norm(prot);
ttx  = ttx/norm(ttx);
tty  = tty/norm(tty);
trot = trot/norm(trot);
tx   = tx/norm(tx);
ty   = ty/norm(ty);
rot  = rot/norm(rot);
for d=1:dim
    pts(idd(d:dim:end),d) = 1;
    tts(itt(d:dim:end),d) = 1;
    pts(:,d)  = pts(:,d)/norm(pts(:,d));
    tts(:,d)  = tts(:,d)/norm(tts(:,d));
end
pert = pert-pts*pts'*pert-prot*prot'*pert;
targ = targ-tts*tts'*targ-trot*trot'*targ;

% on perturbing sites
A = [pts(idd,:)';prot(idd,:)'];
naa = null(A);
Pmat = [naa';A];
dP0 = size(naa,2);
dP1 = size(A,1);
% on target sites
A = [tts(itt,:)';trot(itt,:)'];
naa = null(A);
Tmat = [naa';A];
dT0 = size(naa,2);
dT1 = size(A,1);
% on perturbing and target sites
PTmat = zeros(size(Pmat)+size(Tmat));
PTmat(1:np*dim,1:np*dim) = Pmat;
PTmat(np*dim+(1:nt*dim),np*dim+(1:nt*dim)) = Tmat;
%% Compute Mmatrix
rbb  = zeros(nnb,nnb); % boundary boundary distance
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
    %if dely(n)>ybound/2  % not periodic
    %    dely(n) = dely(n) - ybound;
    %elseif dely(n)<-ybound/2
    %    dely(n) = dely(n) + ybound;
    %end
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
for n1=1:nnb
    for n2=n1:nnb
        dx = bpx(n2)-bpx(n1);
        if dx<-xbound/2
            dx = dx+xbound;
        elseif dx>xbound/2
            dx = dx-xbound;
        end
        dy = bpy(n2)-bpy(n1);
        dr = sqrt(dx^2+dy^2);
        rbb(n1,n2) = dr;
        rbb(n2,n1) = dr;
    end
end
rmin = min(min(rbb(rbb>0)));
rmax = max(max(rbb));
binn = floor((rmax-rmin)/drbin);
rbd  = rmin-drbin/2:drbin:rmin+(binn-1.5)*drbin;
if length(rbd)~=binn
    binn = length(rbd);
end
dn   = max(1,round(2.^(1:log2(nnb/2)/topnb:log2(nnb))-2.^(1-log2(nnb/2)/topnb:log2(nnb/2)/topnb:log2(nnb)-log2(nnb/2)/topnb)));
nbd  = cumsum(dn);
nbd  = nbd(2:end-1);
if topnb~=length(nbd)
    topnb = length(nbd);
end
%Tmatrix = Smatrix';
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
    %if delyw(n)>ybound/2
    %    delyw(n) = delyw(n) - ybound;
    %elseif delyw(n)<-ybound/2
    %    delyw(n) = delyw(n) + ybound;
    %end
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

nl = ceil(nnb/lmem);
pcstt = zeros(1,nl*lmem);
nsp    = round((nnb+3*num)*spf/2);

%% data structure
s = rng('shuffle');
nz     = zeros(temprange,dim*num);
ploc   = zeros(temprange,num);    % local participation ratio
zloc   = zeros(temprange,num);    % local coordination number
bfact  = zeros(temprange,num);    % b-factor
sfact  = zeros(temprange,nnb);    % s-factor
shons  = zeros(temprange,num);    % shear energy on sites
shonb  = zeros(temprange,nnb);    % shear energy on bonds
bkons  = zeros(temprange,num);    % bulk energy on sites
bkonb  = zeros(temprange,nnb);    % bulk energy on bonds
shspr  = zeros(temprange,nnb);    % shear energy on springs
bkspr  = zeros(temprange,nnb);    % bulk energy on springs
cbnr   = zeros(topnb,binn);       % correlation on the number of top sheared bonds and distance r
csnr   = zeros(topnb,binn);       % same correlation but on the springs with finite shear.
nc     = zeros(temprange,1);      % count the number of contributions
cdata  = zeros(temprange,nrec);   % cost
pdata  = zeros(temprange,nrec);   % partition ratio
zdata  = zeros(temprange,nrec);   % coordination of the most participated
estim  = zeros(temprange,nrec);   % energy cost due to stimulus
etarg  = zeros(temprange,nrec);   % energy cost due to target
esttg  = zeros(temprange,nrec);   % energy cost due to cooperativity of stimulus and target
config = zeros(nrec,nnb); %uint32(zeros(nrec,nl));
dd     = zeros(20,nrec);
qq     = zeros(20,nrec);
dd1    = zeros(20,nrec);
qq1    = zeros(20,nrec);
for tt = 1:temprange
    nf(tt) = (min([nflip(tt),nsp,nnb-nsp]));
end

%% Initialize the configurations
cost   = zeros(temprange,1);
prat   = zeros(temprange,1);
coord  = zeros(temprange,1);
occupy = zeros(temprange,nsp-cf);
vacant = zeros(temprange,nnb-nsp);
state  = zeros(temprange,nnb);
%Mmat   = zeros(dim*num);
Qmat   = zeros(dim*num);
for tt = 1:temprange
    occupy(tt,:) = randsample(setdiff(1:nnb,fixs), nsp-cf); %
    vacant(tt,:) = setxor(setdiff(1:nnb,fixs), occupy(tt,:)); %
    state(tt,union(occupy(tt,:),fixs)) = 1+0.5*sign(rand(1,nsp)-0.5); %
    Mmat = kweak*MmatW+Smatrix'*diag(state(tt,:))*Smatrix;
    for i=idd
        Qmat(i,i) = 1;
    end
    Qmat(:,nid) = -Mmat(:,nid);
    fext = Mmat*pert;
    disp = linsolve(Qmat,fext);  %Qmat\fext; %
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext;
    end
    diff = disp(itt)-targ(itt);
    cost(tt) = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
    %compute the coordination for each node
    [cc,~]=hist(reshape(bondnb(union(fixs,occupy(tt,:)),:),1,[]),1:num);
    for d = 1:dim
        nz(tt,d:dim:end) = cc;
    end
    %
    displ = zeros(dim*num,1);
    xid = mod(nid,2)==1;
    yid = mod(nid,2)==0;
    dx = mean(disp(nid(xid)));
    dy = mean(disp(nid(yid)));
    displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
    if nsp>2*xl*yl %subtract the first phonon
        y0 = mean(posy);
        nm = norm(posy-y0);
        dx1 = (posy-y0)*displ(1:2:end)/nm;
        displ(1:2:end) = displ(1:2:end)-dx1*(posy'-y0)/nm;
    end
    prat(tt) = sum(displ(nz(tt,:)>2).^2)^2/num/dim/sum(displ(nz(tt,:)>2).^4);
    coord(tt)= nz(tt,nz(tt,:)>2)*(displ(nz(tt,:)>2).^2)/sum(displ(nz(tt,:)>2).^2);
    %}
end
%% iteration of Monte Carlo
tic;
for t = 1:npara
    for tt = 1:temprange
        %ptarg = targ;
        pbeta = beta(tt);
        pnf   = randi(nf(tt));
        pocc  = occupy(tt,:);
        pvac  = vacant(tt,:);
        pstat = state(tt,:);
        pcost = cost(tt);
        pr    = prat(tt);
        pz    = coord(tt);
        pcdat = zeros(1,prec);
        partr = zeros(1,prec);
        pcoord= zeros(1,prec);
        es    = zeros(1,prec);
        et    = zeros(1,prec);
        est   = zeros(1,prec);
        %partr(1) = prat(tt);
        pconf = zeros(prec,nnb); %uint32(zeros(prec,nl));
        for ii = 1:paraiter

            ctemp = pcost;
            statetemp = pstat;
            % new state
            old = randsample(pocc,pnf);
            if length(pvac)>1
                new = randsample(pvac,pnf);
            else
                new = pvac;
            end
            pstat(new) = pstat(old);
            pstat(old) = 0;
            Mmat = kweak*MmatW+Smatrix'*diag(pstat)*Smatrix;
            % displacement with stimulus
            if allos==1
                Qmat   = zeros(dim*num);
                for i=idd
                    Qmat(i,i) = 1;
                end
                Qmat(:,nid) = -Mmat(:,nid);
                fext = Mmat*pert;
                disp = linsolve(Qmat,fext);  %kweak>0
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext;
                end
                diff = disp(itt)-targ(itt);
                pcost = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
            else
                Qmat = zeros(dim*num);
                for i=idd(1:dP0)
                    Qmat(i,i) = 1;
                end
                Qmat(nid,nid) = -Mmat(nid,nid);
                Qmat(idd,nid) = -Pmat*Mmat(idd,nid);
                Qmat(nid,idd(dP0+1:np*dim)) = -Mmat(nid,idd)*Pmat(dP0+1:np*dim,:)';
                Qmat(idd,idd(dP0+1:np*dim)) = -Pmat*Mmat(idd,idd)*Pmat(dP0+1:np*dim,:)';
                fext = Mmat*pert;
                fext(idd) = Pmat*fext(idd);
                disp = linsolve(Qmat,fext);  %kweak>0
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext; %
                end
                disp(idd(1:dP0)) = 0; % set force to zero
                disp(idd) = Pmat'*disp(idd);
                disp = disp+pert;
                pes  = 0.5*disp'*Mmat*disp;

                % energy of only targ
                Qmat = zeros(dim*num);
                for i=itt(1:dT0)
                    Qmat(i,i) = 1;
                end
                Qmat(nit,nit) = -Mmat(nit,nit);
                Qmat(itt,nit) = -Tmat*Mmat(itt,nit);
                Qmat(nit,itt(dT0+1:nt*dim)) = -Mmat(nit,itt)*Tmat(dT0+1:nt*dim,:)';
                Qmat(itt,itt(dT0+1:nt*dim)) = -Tmat*Mmat(itt,itt)*Tmat(dT0+1:nt*dim,:)';
                fext = Mmat*targ;
                fext(itt) = Tmat*fext(itt);
                disp = linsolve(Qmat,fext);  %kweak>0
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext; %
                end
                disp(itt(1:dT0)) = 0;
                disp(itt) = Tmat'*disp(itt);
                disp = disp+targ;
                pet  = 0.5*disp'*Mmat*disp;

                % displacement with both stimulus and target
                Qmat = zeros(dim*num);
                for i=[idd(1:dP0),itt(1:dT0)]
                    Qmat(i,i) = 1;
                end
                Qmat(nidt,nidt) = -Mmat(nidt,nidt);
                Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
                Qmat(nidt,[idd(dP0+1:np*dim),itt(dT0+1:nt*dim)]) = -Mmat(nidt,idt)*PTmat([dP0+1:np*dim,np*dim+dP0+1:(np+nt)*dim],:)';
                Qmat(idt,[idd(dP0+1:np*dim),itt(dT0+1:nt*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([dP0+1:np*dim,np*dim+dP0+1:(np+nt)*dim],:)';
                fext = Mmat*(pert+targ);
                fext(idt) = PTmat*fext(idt);
                disp = linsolve(Qmat,fext);  %kweak>0
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext; %
                end
                disp([idd(1:dP0),itt(1:dT0)]) = 0;
                disp(idt) = PTmat'*disp(idt);
                disp = disp+pert+targ;
                pest = 0.5*disp'*Mmat*disp;
                pcost = pest-(pes+pet); %cost
            end
            % displacement with translation removed
            dx = mean(disp(1:2:end));
            dy = mean(disp(2:2:end));
            displ = disp-dx*(mod(1:dim*num,2)==1)'-dy*(mod(1:dim*num,2)==0)';
            % metropolis criterion
            if rand(1)<exp(-pbeta*(pcost-ctemp))
                pocc(ismember(pocc,old)) = new;
                pvac(ismember(pvac,new)) = old;
                pr  = sum(displ(nz(tt,:)>2).^2)^2/num/dim/sum(displ(nz(tt,:)>2).^4);
                remv = reshape(bondnb(old,:),1,[]);
                badd = reshape(bondnb(new,:),1,[]);
                for d=1:dim
                    for i=remv
                        nz(tt,dim*(i-1)+d) = nz(tt,dim*(i-1)+d)-1;
                    end
                    for i=badd
                        nz(tt,dim*(i-1)+d) = nz(tt,dim*(i-1)+d)+1;
                    end
                end
                pz  = nz(tt,nz(tt,:)>2)*(displ(nz(tt,:)>2).^2)/sum(displ(nz(tt,:)>2).^2);
            else
                pcost = ctemp;
                pstat = statetemp;
            end
            % record data
            if mod(ii,interv)==0
                pcdat(ii/interv)   = pcost;
                pcstt(1:nnb) = pstat;
                pconf(ii/interv,:) = pstat; %uint32(2.^((lmem-1):-1:0)*reshape(pcstt,lmem,nl));
                partr(ii/interv)   = pr;
                pcoord(ii/interv)  = pz;
                if allos==0
                es(ii/interv)      = pes;
                et(ii/interv)      = pet;
                est(ii/interv)     = pest;
                end
                % compute the eigenvalues
                %{
                pMmat = kweak*MmatW+Smatrix'*diag(pstat)*Smatrix;
                Mmat  = Smatrix'*diag(pstat)*Smatrix;
                [V,D] = eig(full(pMmat));
                [Vr,Dr] = eig(full(Mmat));
                d = diag(D);
                d1= diag(Dr);
                dd(:,ii/interv) = sqrt(d(i0+(1:20)));
                qq(:,ii/interv) = abs(V(:,i0+(1:20))'*displ)/norm(displ);
                dd1(:,ii/interv) = sqrt(d1(i0+(1:20)));
                qq1(:,ii/interv) = abs(Vr(:,i0+(1:20))'*displ)/norm(displ);
                %}
                if t>npara/2
                    nc(tt) = nc(tt)+1;
                    zloc(tt,:) = zloc(tt,:)+nz(tt,1:2:end);
                    ddp = displ(nz(tt,:)>2)';
                    ploc(tt,nz(tt,1:2:end)>2) = ploc(tt,nz(tt,1:2:end)>2)+(ddp(1:2:end).^2+ddp(2:2:end).^2)/sum(ddp.^2);
                    % b-factor and s-factor
                    %{
                    iMmat = pinv(pMmat);
                    bu  = diag(iMmat)';
                    bfact(tt,:) = bfact(tt,:)+bu(1:2:end)+bu(2:2:end);
                    for b = 1:nnb
                        bij = bondnb(b,:);
                        sfact(tt,b) = sfact(tt,b)+bu(2*bij(1)-1)+bu(2*bij(1))+bu(2*bij(2)-1)+bu(2*bij(2))...
                            -iMmat(2*bij(1)-1,2*bij(2)-1)-iMmat(2*bij(2)-1,2*bij(1)-1)-iMmat(2*bij(1),2*bij(2))-iMmat(2*bij(2),2*bij(1));
                    end
                    %}
                    % compute stress tensor
                    %
                    dxx = zeros(1,nnb);
                    dyy = zeros(1,nnb);
                    for b = 1:nnb
                        dxx(b) = displ(2*bondnb(b,1)-1)-displ(2*bondnb(b,2)-1);
                        dyy(b) = displ(2*bondnb(b,1))-displ(2*bondnb(b,2));
                    end
                    exx = delx./delr.^2.*dxx;
                    eyy = dely./delr.^2.*dyy;
                    exy = 0.5*(delx./delr.^2.*dyy+dely./delr.^2.*dxx);
                    etr = 0.5*(exx+eyy);
                    gxx = exx-etr;
                    gyy = eyy-etr;
                    bkonb(tt,:) = exx.^2+eyy.^2;
                    shonb(tt,:) = gxx.^2+gyy.^2+2*exy.^2;
                    bkspr(tt,:) = pstat.*bkonb(tt,:);
                    shspr(tt,:) = pstat.*shonb(tt,:);

                    % measure the two point correlation of the most sheared
                    % links or springs
                    [~,sbord] = sort(shonb(tt,:),'descend'); % order in descend direction
                    [~,ssord] = sort(shspr(tt,:),'descend');
                    for nn = 1:topnb
                        bid = sbord(1:nbd(nn));
                        distr = reshape(rbb(bid,bid),1,[]);
                        rcount= histc(distr(distr>0),rbd);
                        cbnr(nn,:) = cbnr(nn,:)+cumsum(rcount)/nbd(nn)/(nbd(nn)-1);
                        sid = ssord(1:nbd(nn));
                        distr = reshape(rbb(sid,sid),1,[]);
                        rcount= histc(distr(distr>0),rbd);
                        csnr(nn,:) = csnr(nn,:)+cumsum(rcount)/nbd(nn)/(nbd(nn)-1);
                    end
                    %}

                end
            end
        end
        occupy(tt,:) = pocc;
        vacant(tt,:) = pvac;
        state(tt,:)  = pstat;
        cost(tt)     = pcost;
        prat(tt)     = pr;
        coord(tt)    = pz;
        cdata(tt,(t-1)*prec+(1:prec)) = pcdat;
        config((t-1)*prec+(1:prec),:) = pconf;
        pdata(tt,(t-1)*prec+(1:prec)) = partr;
        zdata(tt,(t-1)*prec+(1:prec)) = pcoord;
        estim(tt,(t-1)*prec+(1:prec)) = es;
        etarg(tt,(t-1)*prec+(1:prec)) = et;
        esttg(tt,(t-1)*prec+(1:prec)) = est;
    end
    % switch paratemp data
    if temprange>1
        for tt = 1:temprange-1
            for t1 = 1:temprange-tt
                if rand(1)<exp((beta(t1)-beta(t1+1))*(cost(t1)-cost(t1+1)))
                    occtemp = occupy(t1,:);
                    vactemp = vacant(t1,:);
                    statemp = state(t1,:);
                    costemp = cost(t1);
                    prtemp  = prat(t1);
                    nztemp  = nz(t1,:);
                    ztemp   = coord(t1);
                    occupy(t1,:) = occupy(t1+1,:);
                    vacant(t1,:) = vacant(t1+1,:);
                    state(t1,:)  = state(t1+1,:);
                    cost(t1)     = cost(t1+1);
                    prat(t1)     = prat(t1+1);
                    nz(t1,:)     = nz(t1+1,:);
                    coord(t1)    = coord(t1+1);
                    occupy(t1+1,:) = occtemp;
                    vacant(t1+1,:) = vactemp;
                    state(t1+1,:)  = statemp;
                    cost(t1+1)     = costemp;
                    prat(t1+1)     = prtemp;
                    nz(t1+1,:)     = nztemp;
                    coord(t1+1)    = ztemp;
                end
            end
        end
    end
end
toc;
for tt=1:temprange
    zloc(tt,:) = zloc(tt,:)/nc(tt);
    ploc(tt,:) = ploc(tt,:)/nc(tt);
    bfact(tt,:) = bfact(tt,:)/nc(tt);
    sfact(tt,:) = sfact(tt,:)/nc(tt);
    shons(tt,:) = shons(tt,:)/nc(tt);
    shonb(tt,:) = shonb(tt,:)/nc(tt);
    bkons(tt,:) = bkons(tt,:)/nc(tt);
    bkonb(tt,:) = bkonb(tt,:)/nc(tt);
    shspr(tt,:) = shspr(tt,:)/nc(tt);
    bkspr(tt,:) = bkspr(tt,:)/nc(tt);
end
cbnr = cbnr/nc(1);
csnr = csnr/nc(1);
dirc  = './';
cname = 'Costs';
sname = 'Config';
pname = 'PartRatio';
zname = 'MeanZ';
ppnam = 'PrLocal';
zznam = 'Zlocal';
bfnam = 'Bfactor';
sfnam = 'Sfactor';
ssnam = 'ShearSpring';
sbnam = 'ShearBond';
bsnam = 'BulkSpring';
bbnam = 'BulkBond';
cbnam = 'CorrelBond';
csnam = 'CorrelSpring';
esnam = 'Stimulus';
etnam = 'Target';
stnam = 'Cooperativity';
xname = sprintf('%d', xl);
yname = sprintf('%d', yl);
nsname = sprintf('%d',nsp);
wfname = sprintf('%.4f', kweak);
bname  = sprintf('%.4f',beta_ini);
exname = sprintf('%03d', replica);
costname = [dirc,cname '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(costname, cdata, '\t');
prname = [dirc,pname '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(prname, pdata, '\t');
mzname = [dirc,zname '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(mzname, zdata, '\t');
esname = [dirc,esnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(esname, estim, '\t');
etname = [dirc,etnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(etname, etarg, '\t');
stname = [dirc,stnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(stname, esttg, '\t');
plname = [dirc,ppnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(plname, ploc, '\t');
zlname = [dirc,zznam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(zlname, zloc, '\t');
bfname = [dirc,bfnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(bfname, bfact, '\t');
sfname = [dirc,sfnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(sfname, sfact, '\t');
ssname = [dirc,ssnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(ssname, shspr, '\t');
sbname = [dirc,sbnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(sbname, shonb, '\t');
bsname = [dirc,bsnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(bsname, bkspr, '\t');
bbname = [dirc,bbnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(bbname, bkonb, '\t');
cbname = [dirc,cbnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(cbname,[rbd+0.5*drbin;cbnr], '\t');
csname = [dirc,csnam '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(csname,[rbd+0.5*drbin;csnr], '\t');
%for tt = 1:temprange
%    ttname = sprintf('%d',tt);
    confname = [dirc,sname '_' xname '_' yname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
    dlmwrite(confname, config); %,'precision','%9d');
%end
