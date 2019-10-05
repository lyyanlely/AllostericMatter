function LRCGfix(num,z,f,x0,r)
% this function compute the left right coupling matrices
ff  = 0.001;
dim = 2;
[pos,nb,fl0,fr0,nf]=cutNetworkfix(num,z,f,r-1,x0);
[Smat,~] = SMfix(pos,nb,nf);
pos(:,1) = pos(:,1)-x0;
pos(pos(:,1)<0,1) = pos(pos(:,1)<0,1)+1;
mob = setdiff(1:num,nf);
options.maxit = 1000000;
options.prtlevel = 0;
options.normtol=1e-10; %*ff/sqrt(num)*16;
%% find floppy modes
Z = null(full(Smat'*Smat)); %floppy modes
fl = setdiff(fl0,nf);
fr = setdiff(fr0,nf);
%% projector onto boundary node space
Il = zeros(length(fl)*2,numel(mob)*dim);  
Ir = zeros(length(fr)*2,numel(mob)*dim);
for i=1:length(fl)
    ml = find(mob==fl(i));
    Il(2*i-1,2*ml-1) = 1;
    Il(2*i,2*ml) = 1;
end
for i=1:length(fr)
    mr = find(mob==fr(i));
    Ir(2*i-1,2*mr-1) = 1;
    Ir(2*i,2*mr) = 1;
end
plr= [Il;Ir];  %projector onto the boundary nodes space
%% orthonormalization
pZ = plr*Z;
phi= orth(pZ); % orthonormal base of floppy modes
%% projection operators
Pp = phi*phi';   % projector P
Qp = eye(size(Pp))-Pp;  % projector Q
ll = zeros(size(Pp));
rr = zeros(size(Pp));
nl = 2*length(fl);
nr = 2*length(fr);
for i=1:nl
    ll(i,i)=1;
end
for i=1:nr
    rr(nl+i,nl+i)=1;
end
%% coupling matrix T
Tlr = pinv(ll*Qp*ll)*Pp*rr;
Tlr = Tlr(1:nl,(1:nr)+nl);
Trl = pinv(rr*Qp*rr)*Pp*ll;
Trl = Trl((1:nr)+nl,1:nl);
%% comupte eigenvalues
[~,Sr,Vr] = svd(Trl);
[~,Sl,Vl] = svd(Tlr);
%d  = svd(Tlr*Trl,'econ');
dr = diag(Sr);
dl = diag(Sl);
idr= find(dr>=1);
idl= find(dl>=1);
dr = dr(idr); %(d)>.1);
dl = dl(idl); %(d)>.1);
%% compute the conjugate gradient solution
xx = zeros(length(idr)+length(idl),num);
amp= zeros(length(idr)+length(idl),num);
cdr= zeros(length(idr),1);
cdl= zeros(length(idl),1);
for i = 1:length(idr)
    %pert = zeros(dim*num,1);
    %pert(dim*fl-1,1) = Vr(1:dim:dim*length(fl),idr(i));
    %pert(dim*fl,1) = Vr(2:dim:dim*length(fl),idr(i));
    rid  = setdiff((1:num)',fl);
    pars = initialpars(num,pos,nb,fl);
    options.x0 = pars.x0;
    pars.nvar  = dim*(num-length(fl));
    pars.fgname= 'networkcg';
    pars.pos   = pars.pos+Vr(:,idr(i))/sqrt(num)*ff;
    [x,~,~,~,~]=nlcg(pars,options);
    dx = reshape(x-pars.x0,2,[]);
    xx(i,1:length(rid)) = pars.x0(1:2:2*length(rid))';
    amp(i,1:length(rid))= sqrt(sum(dx.^2));
    cdr(i) = sqrt(num)*norm(reshape(dx(:,ismember(rid,fr)),1,[]));
end
for i = 1:length(idl)
    %pert = zeros(dim*num,1);
    %pert(dim*fr-1,1) = Vl(1:dim:dim*length(fr),idl(i));
    %pert(dim*fr,1) = Vl(2:dim:dim*length(fr),idl(i));
    rid  = setdiff((1:num)',fr);
    pars = initialpars(num,pos,nb,fr);
    options.x0 = pars.x0;
    pars.nvar  = dim*(num-length(fr));
    pars.fgname= 'networkcg';
    pars.pos   = pars.pos+Vl(:,idl(i))/sqrt(num)*ff;
    [x,~,~,~,~]=nlcg(pars,options);
    dx = reshape(x-pars.x0,2,[]);
    xx(i+length(idr),1:length(rid)) = 1-pars.x0(1:2:2*length(rid))';
    amp(i+length(idr),1:length(rid))= sqrt(sum(dx.^2));
    cdl(i) = sqrt(num)*norm(reshape(dx(:,ismember(rid,fl)),1,[]));
end
%[U,S,V] = svd(Trl);
fname = '../Allostery/cg/singval_fix';
aname = '../Allostery/cg/magnif_fix';
pname = '../Allostery/cg/position_fix';
wname = '../Allostery/cg/amplitude_fix';
nname = sprintf('%d',num);
zname = sprintf('%.3f',z);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
fracn = sprintf('%d',-log10(options.normtol));
dtype = '.dat';
filename = [fname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));
filename = [aname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([cdr;cdl]));
filename = [pname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double(xx));
filename = [wname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double(amp));