function [dr,dl,kk,absu] = LRFourier2(num,dim,z,x0,r)
% this function compute the left right coupling matrices
[pos,nb,fl,fr]=cutNetwork(num,dim,z,r,x0);
[Smat,~] = StructureMat(pos,nb);
%% find floppy modes
Z = null(full(Smat'*Smat),'r'); %floppy modes
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
posf = [pos(fl,:);pos(fr,:)];
[Ur,Sr,Vr] = svd(Trl);
[Ul,Sl,Vl] = svd(Tlr);
%d  = svd(Tlr*Trl,'econ');
dr = diag(Sr);
dl = diag(Sl);
dr = dr(dr>=1); %(d)>.1);
dl = dl(dl>=1); %(d)>.1);
%[U,S,V] = svd(Trl);
pZ = plr*Z;
nr = length(dr);
nl = length(dl);
kk = 0:0.1/sqrt(num):1;
absu = zeros(nr+nl,length(kk));
for vi = 1:(nr+nl)
    if vi<=nr
        vecb = [Vr(:,vi);Ur(:,vi)];
        %S  = Sr;
        %ii = vi;
    else
        vecb = [Ul(:,vi-nr);Vl(:,vi-nr)];
        %S  = Sl;
        %ii = vi-nr;
    end
    %coef = pZ\vecb;
    %fmod = Z*coef;
    ff   = zeros(size(vecb));
    ff(1:2:end) = vecb(1:2:end)-mean(vecb(1:2:end));
    ff(2:2:end) = vecb(2:2:end)-mean(vecb(2:2:end));
    %ff(1:2:end) = fmod(1:2:end).*exp(-pos(:,1)*log(S(ii,ii)))-mean(fmod(1:2:end).*exp(-pos(:,1)*log(S(ii,ii))));
    %ff(2:2:end) = fmod(2:2:end).*exp(-pos(:,1)*log(S(ii,ii)))-mean(fmod(2:2:end).*exp(-pos(:,1)*log(S(ii,ii))));
    ff   = ff/norm(ff);
    u = zeros(length(kk),2);
    for ki=1:length(kk)
        u(ki,:) = disfourier(reshape(ff,2,[])',[0;kk(ki)*sqrt(num)],posf);
    end
    u2 = u(:,1).^2+u(:,2).^2;
    absu(vi,:) = abs(sqrt(u2))';
end
%plot(kk,absu)