function LRCGfixV(num,z,f,x0,r)
% this function compute the left right coupling matrices
dim = 2;
[pos,nb,fl0,fr0,nf]=cutNetworkfixV(num,z,f,r-1,x0);
[Smat,~] = SMfixV(pos,nb,nf);
pos(:,1) = pos(:,1)-x0;
pos(pos(:,1)<0,1) = pos(pos(:,1)<0,1)+1;
mob = setdiff(1:num,nf);
options.maxit = 10000;
options.prtlevel = 0;
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
Qp = mp(eye(size(Pp)))-Pp;  % projector Q
ll = mp(zeros(size(Pp)));
rr = mp(zeros(size(Pp)));
nl = 2*length(fl);
nr = 2*length(fr);
for i=1:nl
    ll(i,i)=mp(1);
end
for i=1:nr
    rr(nl+i,nl+i)=mp(1);
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
nfl  = union(fl,nf);
nfr  = union(fr,nf);
lfid = reshape([ismember(nfl,fl);ismember(nfl,fl)],1,[]);
rfid = reshape([ismember(nfr,fr);ismember(nfr,fr)],1,[]);
for i = 1:length(idr)
    rid  = setdiff((1:num)',nfl);
    pars = initialpars(num,pos,nb,nfl);
    options.x0 = pars.x0;
    pars.nvar  = dim*(num-length(nfl));
    pars.fgname= 'networkcg';
    pars.pos(lfid) = pars.pos(lfid)+Vr(:,idr(i))/sqrt(num);
    [x,~,~,~,~]=nlcg(pars,options);
    dx = reshape(x-pars.x0,2,[]);
    xx(i,1:length(rid)) = pars.x0(1:2:2*length(rid))';
    amp(i,1:length(rid))= sqrt(sum(dx.^2));
    cdr(i) = sqrt(num)*norm(reshape(dx(:,ismember(rid,fr)),1,[]));
end
for i = 1:length(idl)
    rid  = setdiff((1:num)',nfr);
    pars = initialpars(num,pos,nb,nfr);
    options.x0 = pars.x0;
    pars.nvar  = dim*(num-length(nfr));
    pars.fgname= 'networkcg';
    pars.pos(rfid) = pars.pos(rfid)+Vl(:,idl(i))/sqrt(num);
    [x,~,~,~,~]=nlcg(pars,options);
    dx = reshape(x-pars.x0,2,[]);
    xx(i+length(idr),1:length(rid)) = 1-pars.x0(1:2:2*length(rid))';
    amp(i+length(idr),1:length(rid))= sqrt(sum(dx.^2));
    cdl(i) = sqrt(num)*norm(reshape(dx(:,ismember(rid,fl)),1,[]));
end
%[U,S,V] = svd(Trl);
fname = 'singval_fix';
aname = 'magnif_fix';
pname = 'position_fix';
wname = 'amplitude_fix';
nname = sprintf('%d',num);
zname = sprintf('%.3f',z);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
fracn = sprintf('%.2f',f);
dtype = '.dat';
filename = [fname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));
filename = [aname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([cdr;cdl]));
filename = [pname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double(xx));
filename = [wname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double(amp));