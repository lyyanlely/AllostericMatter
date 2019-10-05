function [dl,dr,Tlr,Trl] = LRCoupling(lx,ly,dim,z,x0,r)
% this function compute the left right coupling matrices
[pos,nb,fl,fr]=cutNetwork(lx,ly,dim,z,r,x0);
[Smat,~] = StructureMat(pos,nb);
num = size(pos,1);
%% find floppy modes
Z = null(full(Smat'*Smat)); %floppy modes
v1 = zeros(dim*num,1);      % remove translational
v2 = zeros(dim*num,1);
v1(1:dim:end) = 1/sqrt(num);
v2(2:dim:end) = 1/sqrt(num);
Zt = Z-v1*(v1'*Z)-v2*(v2'*Z);
Zn = orth(Zt);
Zn = Zn(:,1:end-2);
%% projector onto boundary node space
Il = zeros(length(fl)*2,numel(pos));  
Ir = zeros(length(fr)*2,numel(pos));
for i=1:length(fl)
    Il(2*i-1,2*fl(i)-1) = 1;
    Il(2*i,2*fl(i)) = 1;
end
for i=1:length(fr)
    Ir(2*i-1,2*fr(i)-1) = 1;
    Ir(2*i,2*fr(i)) = 1;
end
plr= [Il;Ir];  %projector onto the boundary nodes space
%% orthonormalization
pZ = plr*Zn;
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
dr = svd(Trl,'econ');
dl = svd(Tlr,'econ');
%d  = svd(Tlr*Trl,'econ');
dr = dr(dr>=1); %(d)>.1);
dl = dl(dl>=1); %(d)>.1);