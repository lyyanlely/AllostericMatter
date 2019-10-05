function [pos,res] = steepdec(pos,nb,pert)
% compute response field of a perturbation using steepest decent
h   = 0.1;  % step size (viscosity)
eps = 1e-5; % stop criteria
pid = pert~=0;  % where force can be nonzero
%force = zeros(size(pos));
num = size(pos,1);
dim = size(pos,2);

[~,delr,~] = StructureMat(pos,nb);
l0  = delr;
fmod = 1;
cc = 0;

for i = 1:1000
[Smat,~,~] = StructureMat(pos,nb);
[V,D]  = eig(full(Smat'*Smat));
zid = diag(D)<1e-10;
Z  = V(:,zid);
v1 = zeros(dim*num,1);
v2 = zeros(dim*num,1);
v1(1:dim:end) = 1/sqrt(num);
v2(2:dim:end) = 1/sqrt(num);
Zt = Z-v1*(v1'*Z)-v2*(v2'*Z);
Zn = orth(Zt);
Zn = Zn(:,1:end-2);
%nz = size(Zn,2);
%Vn = [Zn';V(:,~zid)'];
pZ = Zn(pid,:);
coef = pinv(pZ)*pert(pid);
fmod = Zn*coef;
norm(fmod)
force = reshape(fmod/norm(fmod),2,[])';
pos   = pos+0.01*force;

[Smat,delr,~] = StructureMat(pos,nb);
str = delr-l0;
res = str*Smat;  % forces on each node

while any(abs(res)>eps)
    %max(res(pni))
    dp  = -res*pinv(full(Smat'*Smat));
    pos = pos+h*reshape(dp,2,[])';
    % recompute the forces
    [Smat,delr,~] = StructureMat(pos,nb);
    str = delr-l0;
    res = str*Smat;
end

end
%pos = pos+pert;