function LRCouplingV(num,x0,r)
% this function compute the left right coupling matrices
%mp.Digits(dps);
dim = 2;
z   = 4.;
[pos,nb,fl,fr]=cutNetworkV(num,dim,z,r,x0);
[Smat,~] = StructureMatV(pos,nb);
%% find floppy modes
Z = null(full(Smat'*Smat)); %floppy modes
%% projector onto boundary node space
Il = mp(zeros(length(fl)*2,numel(pos)));  
Ir = mp(zeros(length(fr)*2,numel(pos)));
for i=1:length(fl)
    Il(2*i-1,2*fl(i)-1) = mp(1);
    Il(2*i,2*fl(i)) = mp(1);
end
for i=1:length(fr)
    Ir(2*i-1,2*fr(i)-1) = mp(1);
    Ir(2*i,2*fr(i)) = mp(1);
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
dr = svd(Trl,'econ');
dl = svd(Tlr,'econ');
%d  = svd(Tlr*Trl,'econ');
dr = dr((dr)>=1);
dl = dl((dl)>=1);

wname = 'singvalue';
nname = sprintf('%d',num);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
dtype = '.dat';
filename = [wname,'_',nname,'_',xname,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));