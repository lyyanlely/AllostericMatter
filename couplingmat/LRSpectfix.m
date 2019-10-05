function LRSpectfix(num,dim,f,x0,r)
% this function compute the left right coupling matrices
[pos,nb,fl0,fr0,nf]=cutNetworkfix(num,dim,f,r,x0);
[Smat,~] = SMfix(pos,nb,nf);
mob = setdiff(1:num,nf);
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
[~,Sr,~] = svd(Trl);
[~,Sl,~] = svd(Tlr);
%d  = svd(Tlr*Trl,'econ');
dr = diag(Sr);
dl = diag(Sl);
dr = dr(dr>=1);
dl = dl(dl>=1);
% spectrum of the amplification
nz = rank(Z);  % number of zero modes
[V,d] = eig(full(Smat'*Smat));
dd = diag(d);
nw = length(d)-nz;
omeg = zeros(1,nw);
magn = zeros(1,nw);
for i = 1:nw
    omeg(i) = sqrt(dd(nz+i));
    phi = plr*V(:,nz+i);
    phi = phi/norm(phi);
    Pp = phi*phi';   % projector P
    Qp = eye(size(Pp))-Pp;  % projector Q
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
fname = '25/singvalue_fix';
wname = '25/spectrum_fix';
nname = sprintf('%d',num);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
fracn = sprintf('%.2f',f);
dtype = '.dat';
filename = [fname,'_',nname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));
filename = [wname,'_',nname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([omeg;magn]));
%[U,S,V] = svd(Trl);
%{
pZ = plr*Z;
nr = length(dr);
nl = length(dl);
kk = 0:0.1/sqrt(num):1;
absu = zeros(nr+nl,length(kk));
for vi = 1:(nr+nl)
    if vi<=nr
        vecb = [Vr(:,vi);Sr(vi,vi)*Ur(:,vi)];
        S  = Sr;
        ii = vi;
    else
        vecb = [Sl(vi-nr,vi-nr)*Ul(:,vi-nr);Vl(:,vi-nr)];
        S  = Sl;
        ii = vi-nr;
    end
    coef = pZ\vecb;
    fmod = Z*coef;
    ff   = zeros(size(fmod));
    ff(1:2:end) = fmod(1:2:end).*exp(-mod(pos(mob,1)+x0,1)*log(S(ii,ii))); %-mean(fmod(1:2:end).*exp(-mod(pos(mob,1)+x0,1)*log(S(ii,ii))));
    ff(2:2:end) = fmod(2:2:end).*exp(-mod(pos(mob,1)+x0,1)*log(S(ii,ii))); %-mean(fmod(2:2:end).*exp(-mod(pos(mob,1)+x0,1)*log(S(ii,ii))));
    ff   = ff/norm(ff);
    u = zeros(length(kk),2);
    for ki=1:length(kk)
        u(ki,:) = disfourier(reshape(ff,2,[])',[0;kk(ki)*sqrt(num)],pos(mob,:));
    end
    u2 = u(:,1).^2+u(:,2).^2;
    absu(vi,:) = abs(sqrt(u2))';
end
%}
%plot(kk,absu)