function LRSpectfixV(num,f,x0,r)
% this function compute the left right coupling matrices
dim = 2;
[pos,nb,fl0,fr0,nf]=cutNetworkfixV(num,dim,f,r-1,x0);
[Smat,~] = SMfixV(pos,nb,nf);
mob = setdiff(1:num,nf);
%% find floppy modes
Z = null((Smat'*Smat)); %floppy modes
fl = setdiff(fl0,nf);
fr = setdiff(fr0,nf);
%% projector onto boundary node space
Il = mp(zeros(length(fl)*2,numel(mob)*dim));  
Ir = mp(zeros(length(fr)*2,numel(mob)*dim));
for i=1:length(fl)
    ml = find(mob==fl(i));
    Il(2*i-1,2*ml-1) = mp(1);
    Il(2*i,2*ml) = mp(1);
end
for i=1:length(fr)
    mr = find(mob==fr(i));
    Ir(2*i-1,2*mr-1) = mp(1);
    Ir(2*i,2*mr) = mp(1);
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
[~,Sr,~] = svd(Trl);
[~,Sl,~] = svd(Tlr);
%d  = svd(Tlr*Trl,'econ');
dr = diag(Sr);
dl = diag(Sl);
dr = dr(dr>=1);
dl = dl(dl>=1);
% spectrum of the amplification
nz = rank(Z);  % number of zero modes
[V,d] = eig((Smat'*Smat));
dd = diag(d);
nw = length(d)-nz;
omeg = mp(zeros(1,nw));
magn = mp(zeros(1,nw));
for i = 1:nw
    omeg(i) = sqrt(dd(nz+i));
    phi = plr*V(:,nz+i);
    phi = phi/norm(phi);
    Pp = phi*phi';   % projector P
    Qp = mp(eye(size(Pp)))-Pp;  % projector Q
    %% coupling matrix T
    Tlr = pinv(ll*Qp*ll)*Pp*rr;
    Tlr = Tlr(1:nl,(1:nr)+nl);
    Trl = pinv(rr*Qp*rr)*Pp*ll;
    Trl = Trl((1:nr)+nl,1:nl);
    %% comupte eigenvalues
    [~,Sr,~] = svd(Trl);
    [~,Sl,~] = svd(Tlr);
    magn(i) = max(Sr(1),Sl(1));
end
fname = 'singvalue_fix';
wname = 'spectrum_fix';
nname = sprintf('%d',num);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
fracn = sprintf('%.2f',f);
dtype = '.dat';
filename = [fname,'_',nname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));
filename = [wname,'_',nname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([omeg;magn]));