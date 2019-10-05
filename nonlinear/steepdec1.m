function [pos,ee] = steepdec1(pos,nb,pert)
% compute response field of a perturbation using steepest decent
h   = 0.1;  % step size (viscosity)
eps = 1e-6; % stop criteria
pid = find(pert==0);  % where force can be nonzero
pp  = union(ceil(pid/2),[]);
%force = zeros(size(pos));
%num = size(pos,1);
%dim = size(pos,2);
ee  = zeros(1,1000);

[~,delr,~] = StructureMat(pos,nb);
l0  = delr;
pos = pos+reshape(pert,2,[])';
[Smat,delr,~] = StructureMat(pos,nb);
str = delr-l0;
res = str*Smat;  % forces on each node
cc = 0;

while any(abs(res)>eps) && cc<=1e6
    cc = cc+1;
    if mod(cc,1000)==0
        str*str'
        ee(cc/1000) = str*str';
    end
    
    force = reshape(res,2,[])';
    pos(pp,:) = pos(pp,:)-h*force(pp,:);
    % recompute the forces
    [Smat,delr,~] = StructureMat(pos,nb);
    str = delr-l0;
    res = str*Smat;
end

%pos = pos+pert;