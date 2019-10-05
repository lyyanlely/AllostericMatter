%function coop3d_ana(replica, allos, beta_ini)
% 3D allosteric model
replica = 1;
xl    = 6;
yl    = 12;
zl    = 12;
spf   = 0.8;
allos = 0;
%dd    = 1;
%dt    = 1.;
kweak = 0.0001;
beta_ini = 0.0001;


np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
xp    = 0; yp = 0; % 0 at the center
xt    = 0; yt = 0; % 0 at the center

eps   = 0.2;   % distort size
warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl*zl;
dim   = 3;
%% define perturb and target displacements
xp  = xl/2+xp;
yp  = yl/2+yp;
xt  = xl/2+xt;
yt  = yl/2+yt;
id  = [(yp-1)*xl+xp,yp*xl+xp+1,(yp+1)*xl+xp,yp*xl+xp];     % indices of the displaced particles
idd = [3*id(1)-2:3*id(1),3*id(2)-2:3*id(2),3*id(3)-2:3*id(3),3*id(4)-2:3*id(4)];  %
it  = [(yt-1)*xl+xt+1,yt*xl+xt+1,(yt+1)*xl+xt+1,yt*xl+xt]+(zl-1)*yl*xl;  % indices of the target particles
itt = [3*it(1)-2:3*it(1),3*it(2)-2:3*it(2),3*it(3)-2:3*it(3),3*it(4)-2:3*it(4)];
idt = [idd,itt];
nid = setxor(1:dim*num,idd);
nit = setxor(1:dim*num,itt);
nidt= setxor(1:dim*num,[idd,itt]);
%ndd = 2*dim+1:dim*num;
pert = zeros(dim*num,1);  % perturb displacement
targ = zeros(dim*num,1);  % target displacement
pts  = zeros(dim*num,dim);  % translation on perturbing nodes
prot = zeros(dim*num,dim);  % rotation on perturbing
tts  = zeros(dim*num,dim);  % translation on target nodes
trot = zeros(dim*num,dim);  % rotation on target
%% construct the embedded network
% build fcc networks from nodes
pos = zeros(dim,num);
for n=1:num
    nz = ceil(n/xl/yl);
    nxy= n - (nz-1)*xl*yl;
    ny = ceil(nxy/xl);
    nx = nxy - (ny-1)*xl;
    pos(3,n) = (nz-1)*sqrt(2)/2;
    pos(2,n) = (ny-1)*sqrt(2)/2;
    pos(1,n) = (nx-1)*sqrt(2)+mod(ny+nz,2)*sqrt(2)/2;
end
% connect the neighbor sites and next neighbors for weak
bondnb = zeros(12*num,2);  % fcc lattice
nnb    = 0;
bondwk = zeros(12*num,2);  % to all next neareat neighbors
nwk    = 0;
fixs = zeros(1,(np+nt)*2);  % labels of springs to be fixed
cf   = 0;
for i=1:num-1
    for j=i+1:num
        dx = zeros(dim,1);
        for d=1:dim
            dx(d) = pos(d,i)-pos(d,j);
        end
        dr = sqrt(sum(dx.^2));
        if dr<1.1
            nnb = nnb+1;
            bondnb(nnb,:) = [i,j];
            % bonds connecting perturbing sites or target sites
            if ismember(i,id)&&ismember(j,id)
                cf = cf+1;
                fixs(cf) = nnb;
            elseif ismember(i,it)&&ismember(j,it)
                cf = cf+1;
                fixs(cf) = nnb;
            end
        elseif dr<1.8
            nwk = nwk+1;
            bondwk(nwk,:) = [i,j];
        end
    end
end
bondnb = bondnb(1:nnb,:);
bondwk = bondwk(1:nwk,:);
nsp    = round(nnb*spf);

dirc  = '../Allostery/3D/';
sname = 'Config3D';
kname = sprintf('%.4f',kweak);
aname = sprintf('%d',allos);
xname  = sprintf('%d',xl);
yname  = sprintf('%d',yl);
zname  = sprintf('%d',zl);
nname  = sprintf('%d',nsp);
tnam = sprintf('%.4f',beta_ini);
rname = sprintf('%03d',replica);
confname = [dirc,sname '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];

if exist(confname,'file')
    config = dlmread(confname);            
%fixs   = fixs(1:cf);
% distort the lattice
for n = 1:num
    nz = ceil(n/xl/yl);
    nxy= n - (nz-1)*xl*yl;
    ny = ceil(nxy/xl);
    nx = nxy - (ny-1)*xl;
    lbl = mod(nx-1+2*mod(ny-1,2)+3*mod(nz-1,2)+2*mod(nx-1,2)*mod(nz-1,2),4);
    if lbl==0
        pos(2,n) = pos(2,n)-eps;
    elseif lbl==1
        pos(1,n) = pos(1,n)+eps;
    elseif lbl==2
        pos(1,n) = pos(1,n)-eps/sqrt(6);
        pos(2,n) = pos(2,n)-eps/sqrt(2);
        pos(3,n) = pos(3,n)+eps/sqrt(3);
    elseif lbl==3
        pos(1,n) = pos(1,n)+eps/sqrt(6);
        pos(2,n) = pos(2,n)+eps/sqrt(2);
        pos(3,n) = pos(3,n)-eps/sqrt(3);
    end
end
%% Compute structure matrix
bp  = zeros(dim,nnb);
del = zeros(dim,nnb);
delr = zeros(1,nnb);
idx = zeros(1,2*dim*nnb);
jdx = zeros(1,2*dim*nnb);
val = zeros(1,2*dim*nnb);
for n = 1:nnb
    ii = bondnb(n,1);
    jj = bondnb(n,2);
    for d = 1:dim
        del(d,n) = pos(d,ii) - pos(d,jj);
        bp(d,n)  = pos(d,jj) + 0.5*del(d,n);
    end
    delr(n) = sqrt(sum(del(:,n).^2));
    for d = 1:dim
        idx(2*dim*(n-1)+d) = n;
        jdx(2*dim*(n-1)+d) = dim*(ii-1)+d;
        val(2*dim*(n-1)+d) = del(d,n)/delr(n);
        idx(2*dim*(n-1)+dim+d) = n;
        jdx(2*dim*(n-1)+dim+d) = dim*(jj-1)+d;
        val(2*dim*(n-1)+dim+d) = -del(d,n)/delr(n);
    end
end
Smatrix = sparse(idx,jdx,val,nnb,dim*num);
% of the weak network
bpw  = zeros(dim,nwk);
delw = zeros(dim,nwk);
delrw = zeros(1,nwk);
idxw = zeros(1,2*dim*nwk);
jdxw = zeros(1,2*dim*nwk);
valw = zeros(1,2*dim*nwk);
for n = 1:nwk
    ii = bondwk(n,1);
    jj = bondwk(n,2);
    for d = 1:dim
        delw(d,n) = pos(d,ii) - pos(d,jj);
        bpw(d,n)  = pos(d,ii)+delw(d,n)/2;
    end
    delrw(n) = sqrt(sum(delw(:,n).^2));
    for d = 1:dim
        idxw(2*dim*(n-1)+d) = n;
        jdxw(2*dim*(n-1)+d) = dim*(ii-1)+d;
        valw(2*dim*(n-1)+d) = delw(d,n)/delrw(n);
        idxw(2*dim*(n-1)+dim+d) = n;
        jdxw(2*dim*(n-1)+dim+d) = dim*(jj-1)+d;
        valw(2*dim*(n-1)+dim+d) = -delw(d,n)/delrw(n);
    end
end
SMatW = sparse(idxw,jdxw,valw,nwk,dim*num);
MmatW = SMatW'*SMatW;

%% initialize the perturbing and target displacements
% perpendicular
signp = -1;
signt = 1;
for i = 1:np
    n0  = id(i);
    nb1 = id(mod(i,np)+1);
    nb2 = id(mod(i+np-2,np)+1);
    v1  = zeros(dim,1);
    v2  = zeros(dim,1);
    for d = 1:dim
        v1(d)  = pos(d,nb1)-pos(d,n0);
        v2(d)  = pos(d,nb2)-pos(d,n0);
    end
    nm  = cross(v1,v2);
    nm  = nm/norm(nm);
    pert(idd((1:dim)+dim*(i-1))) = signp*nm;
    signp = -signp;
end
for i = 1:nt
    n0  = it(i);
    nb1 = it(mod(i,nt)+1);
    nb2 = it(mod(i+nt-2,nt)+1);
    v1  = zeros(dim,1);
    v2  = zeros(dim,1);
    for d = 1:dim
        v1(d)  = pos(d,nb1)-pos(d,n0);
        v2(d)  = pos(d,nb2)-pos(d,n0);
    end
    nm  = cross(v1,v2);
    nm  = nm/norm(nm);
    targ(itt((1:dim)+dim*(i-1))) = signt*nm;
    signt = -signt;
end
%}
% translational and rotational parts at the stimulus site and target site
pcen = [mean(pos(1,id)),mean(pos(2,id)),mean(pos(3,id))]; % center of perturbing nodes
tcen = [mean(pos(1,it)),mean(pos(2,it)),mean(pos(3,it))]; % center of targeting nodes
for d=1:dim
    pts(idd(d:dim:end),d) = 1;
    tts(itt(d:dim:end),d) = 1;
    for i = 1:np
        dx = pos(mod(d,dim)+1,id(i)) - pcen(mod(d,dim)+1);
        dy = pos(mod(d+1,dim)+1,id(i)) - pcen(mod(d+1,dim)+1);
        prot(dim*id(i)-dim+mod(d,dim)+1,d)  = -dy;
        prot(dim*id(i)-dim+mod(d+1,dim)+1,d)= dx;
    end
    for i = 1:nt
        dx = pos(mod(d,dim)+1,it(i)) - tcen(mod(d,dim)+1);
        dy = pos(mod(d+1,dim)+1,it(i)) - tcen(mod(d+1,dim)+1);
        trot(dim*it(i)-dim+mod(d,dim)+1,d)  = -dy;
        trot(dim*it(i)-dim+mod(d+1,dim)+1,d)= dx;
    end
    pts(:,d)  = pts(:,d)/norm(pts(:,d));
    prot(:,d) = prot(:,d)/norm(prot(:,d));
    tts(:,d)  = tts(:,d)/norm(tts(:,d));
    trot(:,d) = trot(:,d)/norm(trot(:,d));
end
prot = orth(prot);
trot = orth(trot);
pert = pert-pts*pts'*pert-prot*prot'*pert;
targ = targ-tts*tts'*targ-trot*trot'*targ;

% construct the transform matrix to consider the translation and rotations
% on perturbing sites
A = [pts(idd,:)';prot(idd,:)'];
naa = null(A);
Pmat = [naa';A];
% on target sites
A = [tts(itt,:)';trot(itt,:)'];
naa = null(A);
Tmat = [naa';A];
% on perturbing and target sites
PTmat = zeros(size(Pmat)+size(Tmat));
PTmat(1:np*dim,1:np*dim) = Pmat;
PTmat(np*dim+(1:nt*dim),np*dim+(1:nt*dim)) = Tmat;
%% data structure

%minp = xl/num;
neig = ceil(num*dim/20);
nconf = size(config,1);
spf  = nsp/nnb;
nend = nconf;  % number of configurations each sequence
ns  = 1; %nconf/nend;  % number of sequences
eps = 1e-10;
%{
sid = [];
for i=0:(xl-np)
    sid = [sid;(1:np)+i];
end
sid = [sid;it];
nsid = size(sid,1);

%pos = [posx;posy]
%}


s = rng('shuffle');
msig  = zeros(ns,nnb);  % mean occupancy
flip  = zeros(ns,nnb);  % cost of flipping
mdisp = zeros(ns,num*dim);
magn  = zeros(ns,num);
cons  = zeros(ns,num);  % conservation on each node
fcst  = zeros(ns,num);  % flipping cost on each node
eshr  = zeros(ns,num);  % shear pseudoenergy
sbfc  = zeros(ns,num);  % strain B-factor
%sencost = zeros(ns,nsid);
%senmagn = zeros(ns,nsid*num);
dd = zeros(ns,neig*nend);
pr = zeros(ns,neig*nend);
qq = zeros(ns,neig*nend);
rd = zeros(ns,neig*nend);  % random configurations
rp = zeros(ns,neig*nend);
rq = zeros(ns,neig*nend);
for r=1:ns
    tic;
    %% conservation
    msig(r,:) = mean(config((r-1)*nend+1:r*nend,:));
    msigc   = min(max(msig(r,:),eps),1-eps);
    conserv = msigc.*log(msigc/spf)+(1-msigc).*log((1-msigc)/(1-spf));  % conservation on each link
    for i=1:nend
        state = config((r-1)*nend+i,:);
        %% response
        Mmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
        disp = computeResp(Mmat,pos,allos,pert,idd,Pmat);
        nd   = norm(disp);
        disp = disp'/nd;
        mdisp(r,:) = mdisp(r,:)+disp/nend;
        mm = zeros(1,num);
        for d = 1:dim
            mm = mm+disp(d:dim:end).^2;
        end
        magn(r,:) = magn(r,:)+sqrt(mm)/nend*sqrt(num);
        %% pseudo-energy
        [strn,vor] = computeStrain(disp',pos,bondnb);
        [bk,sh] = strainEnergy(strn);
        eshr(r,:) = eshr(r,:)+sh'/nend;
        %% omega, participation, overlap
        [V,d]=eigs(Mmat+eps*eye(num*dim),neig,'sm');
        pr(r,(i-1)*neig+1:i*neig) = computePart(V,2);
        dd(r,(i-1)*neig+1:i*neig) = diag(d)';
        qq(r,(i-1)*neig+1:i*neig) = (disp*V).^2;
        %% B-factor
        if i==nend
            %tic;
            [bf,stb]  = compBfactor(Mmat,pos,bondnb);
            %toc;
            sbfc(r,:) = stb';
        end
        %% sensitivity test
        %{
        if i==nend
            %tic;
            for ii=1:nsid
                [pp,tt,ids,mats] = prepareTest(pos,sid(ii,:),it);
                [pcost,displ] = comppcost(Mmat,pp,tt,allos,pos,ids,mats);
                displ = displ'/nd;
                sencost(r,ii) = pcost;
                mm = zeros(1,num);
                for d = 1:dim
                    mm = mm+displ(d:dim:end).^2;
                end
                senmagn(r,(ii-1)*num+1:ii*num) = sqrt(mm)*sqrt(num);
            end
            %toc;
        end
        %}
        %% flipping cost
        if i==nend
            %tic;
            [pcost,~] = comppcost(Mmat,pert,targ,allos,pos,{idd,nid,itt,nit,idt,nidt},{Pmat,Tmat,PTmat});
            pcost0 = pcost;
            for n=1:nnb
                state1 = state;
                state1(n) = mod(state(n)+1,2);
                Mmat = kweak*MmatW+Smatrix'*diag(state1)*Smatrix;
                [pcost,~] = comppcost(Mmat,pert,targ,allos,pos,{idd,nid,itt,nit,idt,nidt},{Pmat,Tmat,PTmat});
                flip(r,n) = pcost0-pcost;
            end
            %toc;
        end
        %% random configuration
        occ = zeros(nnb,1);
        idx = randsample(nnb,nsp);
        occ(idx) = 1;
        Mmat = kweak*MmatW+Smatrix'*diag(occ)*Smatrix;
        [V,d]=eigs(Mmat+eps*eye(num*dim),neig,'sm');
        rp(r,(i-1)*neig+1:i*neig) = computePart(V,2);
        rd(r,(i-1)*neig+1:i*neig) = diag(d)';
        rq(r,(i-1)*neig+1:i*neig) = (disp*V).^2;
    end
    % from bond to site
    for n=1:num   % average link values to each node
        flg = ismember(bondnb(:,1),n)|ismember(bondnb(:,2),n);
        cons(r,n) = mean(conserv(flg));
        fcst(r,n) = mean(flip(r,flg));
    end
    toc;
end

dpnam = 'MDisp3D';
bevon = 'Evobond3D';
sfeat = 'SiteFeat3D';
specn = 'Spectrum3D';
mdpname = [dpnam '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
bevname = [bevon '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
sftname = [sfeat '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
spcname = [specn '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
dlmwrite(mdpname,mdisp,'\t');
dlmwrite(bevname,[msig;flip],'\t');
dlmwrite(sftname,[magn;cons;fcst;eshr;sbfc],'\t');
dlmwrite(spcname,[dd;pr;qq;rd;rp;rq],'\t');
end