%function allostery3(replica, xl, yl, zl, spf, beta_ini)
% 3D allosteric model
replica = 1;
xl    = 6;
yl    = 12;
zl    = 12;
spf   = 0.83;
allos = 0;
%dd    = 1;
%dt    = 1.;
kweak = 0.0001;
beta_ini = 0.01;

np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
xp    = 0; yp = 0; % 0 at the center
xt    = 0; yt = 0; % 0 at the center

iteration = 100000;
interv = 1000;
nrec   = round(iteration/interv);

eps   = 0.;   % distort size

warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl*zl;
dim   = 3;
beta  = 1/beta_ini;

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
ndd = 2*dim+1:dim*num;
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
fixs   = fixs(1:cf);
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
%pert(dim*id) = [-1,1,-1,1]*dd;
%targ(dim*it) = [1,-1,1,-1]*dt;
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
%A = eye(dim*num);
%A(:,idd(2*dim+(1:dim))) = pts;
%A(:,idd(3*dim+(1:dim))) = prot;
%A(:,idd(1:2*dim)) = [];
%naa= null(A');
%Pmat = [naa';A'];
A = [pts(idd,:)';prot(idd,:)'];
naa = null(A);
Pmat = [naa';A];
% on target sites
A = [tts(itt,:)';trot(itt,:)'];
naa = null(A);
Tmat = [naa';A];
% on perturbing and target sites
%A = [pts(idt,:)';prot(idt,:)';tts(idt,:)';trot(idt,:)'];
%naa = null(A);
PTmat = zeros(size(Pmat)+size(Tmat));
PTmat(1:np*dim,1:np*dim) = Pmat;
PTmat(np*dim+(1:nt*dim),np*dim+(1:nt*dim)) = Tmat;
%A = eye(dim*num);
%A(:,idd(2*dim+(1:dim))) = pts;
%A(:,idd(3*dim+(1:dim))) = prot;
%A(:,itt(2*dim+(1:dim))) = tts;
%A(:,itt(3*dim+(1:dim))) = trot;
%A(:,[idd(1:2*dim),itt(1:2*dim)]) = [];
%naa= null(A');
%PTmat = [naa';A'];
%% data structure
nsp    = round(nnb*spf);
s = rng('shuffle');
cdata  = zeros(1,nrec);   % cost
zdata  = zeros(nrec,num); % coordination number
estim  = zeros(1,nrec);   % energy cost due to stimulus
etarg  = zeros(1,nrec);   % energy cost due to target
esttg  = zeros(1,nrec);   % energy cost due to cooperativity of stimulus and target
config = zeros(nrec/2,nnb);

%% Initialize the configurations
state  = zeros(1,nnb);

occupy = randsample(setdiff(1:nnb,fixs), nsp-cf); %
vacant = setxor(setdiff(1:nnb,fixs), occupy); %
state(union(occupy,fixs)) = 1; %
Mmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
if allos == 1  %allosteric cost
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
    cost = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
else % cooperative cost
    % energy of only pert
    Qmat = zeros(dim*num);
    for i=idd(1:2*dim)
        Qmat(i,i) = 1;
    end
    Qmat(nid,nid) = -Mmat(nid,nid);
    Qmat(idd,nid) = -Pmat*Mmat(idd,nid);
    Qmat(nid,idd(2*dim+1:4*dim)) = -Mmat(nid,idd)*Pmat(2*dim+1:4*dim,:)';
    Qmat(idd,idd(2*dim+1:4*dim)) = -Pmat*Mmat(idd,idd)*Pmat(2*dim+1:4*dim,:)';
    fext = Mmat*pert;
    fext(idd) = Pmat*fext(idd);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(idd(1:2*dim)) = 0; % set force to zero
    disp(idd) = Pmat'*disp(idd);
    disp = disp+pert;
    pes  = 0.5*disp'*Mmat*disp;
    
    % energy of only targ
    Qmat = zeros(dim*num);
    for i=itt(1:2*dim)
        Qmat(i,i) = 1;
    end
    Qmat(nit,nit) = -Mmat(nit,nit);
    Qmat(itt,nit) = -Tmat*Mmat(itt,nit);
    Qmat(nit,itt(2*dim+1:4*dim)) = -Mmat(nit,itt)*Tmat(2*dim+1:4*dim,:)';
    Qmat(itt,itt(2*dim+1:4*dim)) = -Tmat*Mmat(itt,itt)*Tmat(2*dim+1:4*dim,:)';
    fext = Mmat*targ;
    fext(itt) = Tmat*fext(itt);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(itt(1:2*dim)) = 0;
    disp(itt) = Tmat'*disp(itt);
    disp = disp+targ;
    pet  = 0.5*disp'*Mmat*disp;
    %{
    pMmat = Tmat*Mmat*Tmat';
    Qmat(:,ndd) = -pMmat(:,ndd);
    fext = Tmat*(Mmat*targ);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(1:2*dim) = 0;
    disp = Tmat'*disp+targ;
    pet = 0.5*disp'*Mmat*disp;
    %}
    
    % displacement with both stimulus and target
    Qmat = zeros(dim*num);
    for i=[idd(1:2*dim),itt(1:2*dim)]
        Qmat(i,i) = 1;
    end
    Qmat(nidt,nidt) = -Mmat(nidt,nidt);
    Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
    Qmat(nidt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)]) = -Mmat(nidt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
    Qmat(idt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
    fext = Mmat*(pert+targ);
    fext(idt) = PTmat*fext(idt);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp([idd(1:2*dim),itt(1:2*dim)]) = 0;
    disp(idt) = PTmat'*disp(idt);
    disp = disp+pert+targ;
    pest = 0.5*disp'*Mmat*disp;
    %{
    pMmat = PTmat*Mmat*PTmat';
    Qmat(:,4*dim+1:dim*num) = -pMmat(:,4*dim+1:dim*num);
    fext = PTmat*(Mmat*(pert+targ));
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(1:4*dim) = 0;
    disp = PTmat'*disp+pert+targ;
    pest = 0.5*disp'*Mmat*disp;
    %}
    cost = -pes-pet+pest; %cost
end
%compute the coordination for each node
[cc,~] = hist(reshape(bondnb(union(fixs,occupy),:),1,[]),1:num);
ncoord = cc;
%
%% iteration of Monte Carlo
tic;
nf   = randi(1);
es    = zeros(1,nrec);
et    = zeros(1,nrec);
est   = zeros(1,nrec);
for ii = 1:iteration
    
    ctemp = cost;
    statetemp = state;
    % new state
    old = randsample(occupy,nf);
    if length(vacant)>1
        new = randsample(vacant,nf);
    else
        new = vacant;
    end
    state(old) = 0;
    state(new) = 1;
    Mmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
    % compute new cost
    if allos == 1  %allosteric cost
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
        cost = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
    else % cooperative cost
        % energy of only pert
        Qmat = zeros(dim*num);
        for i=idd(1:2*dim)
            Qmat(i,i) = 1;
        end
        Qmat(nid,nid) = -Mmat(nid,nid);
        Qmat(idd,nid) = -Pmat*Mmat(idd,nid);
        Qmat(nid,idd(2*dim+1:4*dim)) = -Mmat(nid,idd)*Pmat(2*dim+1:4*dim,:)';
        Qmat(idd,idd(2*dim+1:4*dim)) = -Pmat*Mmat(idd,idd)*Pmat(2*dim+1:4*dim,:)';
        fext = Mmat*pert;
        fext(idd) = Pmat*fext(idd);
        disp = linsolve(Qmat,fext);  %kweak>0
        if any(isnan(disp)|abs(disp)>1e5)
            disp = pinv(Qmat)*fext; %
        end
        disp(idd(1:2*dim)) = 0; % set force to zero
        disp(idd) = Pmat'*disp(idd);
        disp = disp+pert;
        pes  = 0.5*disp'*Mmat*disp;
        
        % energy of only targ
        Qmat = zeros(dim*num);
        for i=itt(1:2*dim)
            Qmat(i,i) = 1;
        end
        Qmat(nit,nit) = -Mmat(nit,nit);
        Qmat(itt,nit) = -Tmat*Mmat(itt,nit);
        Qmat(nit,itt(2*dim+1:4*dim)) = -Mmat(nit,itt)*Tmat(2*dim+1:4*dim,:)';
        Qmat(itt,itt(2*dim+1:4*dim)) = -Tmat*Mmat(itt,itt)*Tmat(2*dim+1:4*dim,:)';
        fext = Mmat*targ;
        fext(itt) = Tmat*fext(itt);
        disp = linsolve(Qmat,fext);  %kweak>0
        if any(isnan(disp)|abs(disp)>1e5)
            disp = pinv(Qmat)*fext; %
        end
        disp(itt(1:2*dim)) = 0;
        disp(itt) = Tmat'*disp(itt);
        disp = disp+targ;
        pet  = 0.5*disp'*Mmat*disp;
        
        % displacement with both stimulus and target
        Qmat = zeros(dim*num);
        for i=[idd(1:2*dim),itt(1:2*dim)]
            Qmat(i,i) = 1;
        end
        Qmat(nidt,nidt) = -Mmat(nidt,nidt);
        Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
        Qmat(nidt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)]) = -Mmat(nidt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
        Qmat(idt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
        fext = Mmat*(pert+targ);
        fext(idt) = PTmat*fext(idt);
        disp = linsolve(Qmat,fext);  %kweak>0
        if any(isnan(disp)|abs(disp)>1e5)
            disp = pinv(Qmat)*fext; %
        end
        disp([idd(1:2*dim),itt(1:2*dim)]) = 0;
        disp(idt) = PTmat'*disp(idt);
        disp = disp+pert+targ;
        pest = 0.5*disp'*Mmat*disp;
        cost = -pes-pet+pest; %cost
    end
    % metropolis criterion
    if rand(1)<exp(-beta*(cost-ctemp))
        occupy(ismember(occupy,old)) = new;
        vacant(ismember(vacant,new)) = old;
        remv = reshape(bondnb(old,:),1,[]);
        badd = reshape(bondnb(new,:),1,[]);
        for i=remv
            ncoord(i) = ncoord(i)-1;
        end
        for i=badd
            ncoord(i) = ncoord(i)+1;
        end
    else
        cost  = ctemp;
        state = statetemp;
    end
    % record data
    if mod(ii,interv)==0
        cdata(ii/interv)   = cost;
        zdata(ii/interv,:) = ncoord;
        if allos~=1
            estim(ii/interv)   = pes;
            etarg(ii/interv)   = pet;
            esttg(ii/interv)   = pest;
        end
        if ii>iteration/2
            config(ii/interv-nrec/2,1:nnb) = state;
        end
    end
end
toc;
cname = 'Costs3D';
sname = 'Config3D';
cdnam = 'Zlocal3D';
esnam = 'Stimulus3D';
etnam = 'Target3D';
stnam = 'Cooperativity3D';
xname = sprintf('%d', xl);
yname = sprintf('%d', yl);
zname = sprintf('%d', zl);
nsname = sprintf('%d',nsp);
wfname = sprintf('%.4f', kweak);
bname  = sprintf('%.4f',beta_ini);
exname = sprintf('%03d', replica);
costname = [cname '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(costname, cdata, '\t');
mzname = [cdnam '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(mzname, zdata, '\t');
confname = [sname '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(confname, config,'\t');
if allos~=1
    esname = [esnam '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
    dlmwrite(esname, estim, '\t');
    etname = [etnam '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
    dlmwrite(etname, etarg, '\t');
    stname = [stnam '_' xname '_' yname '_' zname '_' nsname '_' bname '_' wfname '_' exname '.dat'];
    dlmwrite(stname, esttg, '\t');
end