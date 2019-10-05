function LRMagnfix(num,z,f,x0,r)
% this function compute the left right coupling matrices
%num = 64;
alpha = 100*eps;
dim = 2;
nn  = round(sqrt(num)); % number of random trials
dr  = zeros(nn,1);
dl  = zeros(nn,1);
%z = 4.125;
%f = 0.0625;
%r = 1;
%x0 = 0.5;
[pos,nb,fl0,fr0,nf]=cutNetworkfix(num,z,f,r-1,x0);
[Smat,~] = SMfix(pos,nb,nf);
mob = setdiff(1:num,nf);
nm  = length(mob);
xx  = zeros(1,2*nn*nm);  % distance to the imposed boundary
amp = zeros(1,2*nn*nm);  % amplitude of displacements
%% find floppy modes
rng('shuffle');
fl = setdiff(fl0,nf);
fr = setdiff(fr0,nf);
%% compute response of a random perturbation
id = []; it = [];
for i = 1:length(fl)
    i0 = find(mob==fl(i));
    id = [id,2*i0-1,2*i0];
end
for i = 1:length(fr)
    i0 = find(mob==fr(i));
    it = [it,2*i0-1,2*i0];
end
nid  = setdiff(1:numel(mob)*dim,id);
Qmat = zeros(numel(mob)*dim);
for i=id
    Qmat(i,i) = 1;
end
Mmat = alpha*eye(numel(mob)*dim)+Smat'*Smat;
Qmat(:,nid) = -Mmat(:,nid);
for n = 1:nn
pert = zeros(numel(mob)*dim,1);
pert(id) = randn(2*length(fl),1);
pert(id(1:2:end)) = pert(id(1:2:end))-mean(pert(id(1:2:end)));
pert(id(2:2:end)) = pert(id(2:2:end))-mean(pert(id(2:2:end)));
pert = pert/norm(pert(id));
fext = Mmat*pert;
disp = linsolve(Qmat,fext);
disp(id) = pert(id);  % displacement into imposed
dr(n) = norm(disp(it));
xx((n-1)*nm+(1:nm))  = rem(pos(mob,1)+1-x0,1); %relative distance to left boundary
amp((n-1)*nm+(1:nm)) = sqrt(disp(1:2:end).^2+disp(2:2:end).^2); % amplitude of displacement
end
%% comupte response on the other side
id = []; it = [];
for i = 1:length(fr)
    i0 = find(mob==fr(i));
    id = [id,2*i0-1,2*i0];
end
for i = 1:length(fl)
    i0 = find(mob==fl(i));
    it = [it,2*i0-1,2*i0];
end
nid  = setdiff(1:numel(mob)*dim,id);
Qmat = zeros(numel(mob)*dim);
for i=id
    Qmat(i,i) = 1;
end
Mmat = alpha*eye(numel(mob)*dim)+Smat'*Smat;
Qmat(:,nid) = -Mmat(:,nid);
for n = 1:nn
pert = zeros(numel(mob)*dim,1);
pert(id) = randn(2*length(fr),1);
pert(id(1:2:end)) = pert(id(1:2:end))-mean(pert(id(1:2:end)));
pert(id(2:2:end)) = pert(id(2:2:end))-mean(pert(id(2:2:end)));
pert = pert/norm(pert(id));
fext = Mmat*pert;
disp = linsolve(Qmat,fext);
disp(id) = pert(id);
dl(n) = norm(disp(it));
xx((nn+n-1)*nm+(1:nm))  = rem(1+x0-pos(mob,1),1); %relative distance to right boundary
amp((nn+n-1)*nm+(1:nm)) = sqrt(disp(1:2:end).^2+disp(2:2:end).^2); % amplitude of displacement
end
%% fourier transform
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
fname = '../Allostery/magn/magnif_fix';
wname = '../Allostery/magn/amplitude_fix';
nname = sprintf('%d',num);
zname = sprintf('%.3f',z);
xname = sprintf('%.2f',x0);
rname = sprintf('%03d',r);
fracn = sprintf('%.2f',f);
dtype = '.dat';
filename = [fname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([dr;dl]));
filename = [wname,'_',nname,'_',zname,'_',xname,'_',fracn,'_',rname,dtype];
dlmwrite(filename,double([xx;amp]));