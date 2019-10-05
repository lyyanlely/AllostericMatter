function displ = rmvtsrot(disp,pos,opt)
if nargin<3
    opt.bd = 1;  % 0 for periodic boundary
end
dsz = size(pos);
dim = dsz(1);
num = dsz(2);
trans = zeros(dim*num,dim);
drot = round(dim*(dim-1)/2);
rotat = zeros(dim*num,drot);
for d = 1:dim
    trans(d:dim:end,d) = 1;
end
for d = 1:dim
    trans(:,d) = trans(:,d)/norm(trans(:,d));
end
cen = zeros(dim,1);
for d=1:dim
    cen(d) = mean(pos(d,:));
end
for d=1:drot
    d1 = mod(d,dim)+1;
    d2 = mod(d+1,dim)+1;
    dx = pos(d1,:)-cen(d1);
    dy = pos(d2,:)-cen(d2);
    rotat(d1:dim:end,d) = -dy;
    rotat(d2:dim:end,d) = dx;
end
for d=1:drot
    rotat(:,d) = rotat(:,d)/norm(rotat(:,d));
end
displ=disp;
for d=1:dim
    displ=displ-trans(:,d)'*displ*trans(:,d);
end
if opt.bd
    for d=1:drot
        displ=displ-rotat(:,d)'*displ*rotat(:,d);
    end
end