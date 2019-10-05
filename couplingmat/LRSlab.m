function [dl,dr,Tlr,Trl] = LRSlab(num,dim,z,xn,r)
% this function compute the left right coupling matrices
[pos,nb,fl,fr]=cutSlab(num,dim,z,xn,r);
np  = length(pos);
Tlr = cell(np,1);
Trl = cell(np,1);
dl  = [];
dr  = [];
for i=1:np
[Smat,~] = StructureMat(pos{i},nb{i});
%% find floppy modes
Z = null(full(Smat'*Smat)); %floppy modes
v1 = zeros(numel(pos{i}),1);      % remove translational
v2 = zeros(numel(pos{i}),1);
v1(1:dim:end) = 1/sqrt(size(pos{i},1));
v2(2:dim:end) = 1/sqrt(size(pos{i},1));
Zt = Z-v1*(v1'*Z)-v2*(v2'*Z);
Zn = orth(Zt);
Zn = Zn(:,1:end-2);
%% projector onto boundary node space
Il = zeros(length(fl{i})*2,numel(pos{i}));  
Ir = zeros(length(fr{i})*2,numel(pos{i}));
for ii=1:length(fl{i})
    Il(2*ii-1,2*fl{i}(ii)-1) = 1;
    Il(2*ii,2*fl{i}(ii)) = 1;
end
for ii=1:length(fr{i})
    Ir(2*ii-1,2*fr{i}(ii)-1) = 1;
    Ir(2*ii,2*fr{i}(ii)) = 1;
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
nl = 2*length(fl{i});
nr = 2*length(fr{i});
for ii=1:nl
    ll(ii,ii)=1;
end
for ii=1:nr
    rr(nl+ii,nl+ii)=1;
end
%% coupling matrix T
Tlr{i} = pinv(ll*Qp*ll)*Pp*rr;
Tlr{i} = Tlr{i}(1:nl,(1:nr)+nl);
Trl{i} = pinv(rr*Qp*rr)*Pp*ll;
Trl{i} = Trl{i}((1:nr)+nl,1:nl);
%% comupte eigenvalues
dri = svd(Trl{i},'econ');
dli = svd(Tlr{i},'econ');
%d  = svd(Tlr{i}*Trl{i});
dr = [dr;dri(dri>1)];
dl = [dl;dli(dli>1)];
end