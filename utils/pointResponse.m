function [disp,pert,fext] = pointResponse(point,Smatrix,pos)
% this function compute the point response to a dilation
num = size(pos,1);
dim = size(pos,2);
ds  = sqrt((pos(:,1)-point(1)).^2+(pos(:,2)-point(2)).^2);
[~,ord] = sort(ds);
id1 = ord(1);
id2 = ord(2);
id  = [2*id1-1,2*id1,2*id2-1,2*id2];
nid = setdiff(1:num*dim,id);
pert = zeros(num*dim,1);

dx  = pos(id2,1)-pos(id1,1);
dy  = pos(id2,2)-pos(id1,2);
while dx>0.5
    dx = dx-1;
end
while dx<-0.5
    dx = dx+1;
end
while dy>0.5
    dy = dy-1;
end
while dy<-0.5
    dy = dy+1;
end
dr  = sqrt(dx^2+dy^2);

pert(id(1)) = -dx/2;
pert(id(2)) = -dy/2;
pert(id(3)) = dx/2;
pert(id(4)) = dy/2;

Mmat = full(Smatrix'*Smatrix);
Qmat = zeros(dim*num);
for i=id
    Qmat(i,i) = 1;
end
Qmat(:,nid) = -Mmat(:,nid);
%fext = zeros(num*dim,1);
%fext(id(1)) = -dx/dr/2;
%fext(id(2)) = -dy/dr/2;
%fext(id(3)) = dx/dr/2;
%fext(id(4)) = dy/dr/2;
fext = Mmat*pert;
disp = linsolve(Qmat,fext);  %
if any(isnan(disp)|abs(disp)>1e5)
    disp = -pinv(Qmat)*fext; %
end
disp(1:2:end) = disp(1:2:end)-mean(disp(1:2:end));
disp(2:2:end) = disp(2:2:end)-mean(disp(2:2:end));