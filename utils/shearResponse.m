function [disp,pert,fext] = shearResponse(upp,low,Smatrix,pos,fg)
% this function compute the point response to a shear
warning off
a = 0.25;
num = size(pos,1);
dim = size(pos,2);
dsu = sqrt((pos(:,1)-upp(1)).^2+(pos(:,2)-upp(2)).^2);
dsl = sqrt((pos(:,1)-low(1)).^2+(pos(:,2)-low(2)).^2);
if nargin<5
    fg=0;
end
if fg==0
    inn = find(pos(:,1)>low(1)-a/sqrt(num)&pos(:,1)<upp(1)+a/sqrt(num)|pos(:,2)<low(2)+a/sqrt(num)|pos(:,2)<upp(2)+a/sqrt(num));
elseif fg==1
    inn = find(pos(:,1)>low(1)-a/sqrt(num)&pos(:,1)<upp(1)+a/sqrt(num)|pos(:,2)>low(2)-a/sqrt(num)|pos(:,2)>upp(2)-a/sqrt(num));
elseif fg==2
    in1 = find(pos(:,1)>low(1)-a/sqrt(num)&pos(:,1)<upp(1)+a/sqrt(num)|pos(:,2)<low(2)+a/sqrt(num)|pos(:,2)<upp(2)+a/sqrt(num));
    in2 = find(pos(:,1)>low(1)-a/sqrt(num)&pos(:,1)<upp(1)+a/sqrt(num)|pos(:,2)>low(2)-a/sqrt(num)|pos(:,2)>upp(2)-a/sqrt(num));
end
[~,odu] = sort(dsu);
[~,odl] = sort(dsl);
if fg<2
    i = 1;
    id1 = odu(i);
    while ismember(id1,inn)
        i = i+1;
        id1 = odu(i);
    end
    i = 1;
    id2 = odl(i);
    while ismember(id2,inn)
        i = i+1;
        id2 = odl(i);
    end
    id  = [2*id1-1,2*id1,2*id2-1,2*id2];
    nid = setdiff(1:num*dim,id);
    ptx = zeros(num*dim,1);
    pty = zeros(num*dim,1);
    pert = zeros(num*dim,1);
    
    pert(id(1)) = -1; pert(id(2)) = 0; pert(id(3)) = 1; pert(id(4)) = 0;
    if fg==1
        pert = -pert;
    end
    ptx(id(1:2:end)) = 1/sqrt(2);
    pty(id(2:2:end)) = 1/sqrt(2);
    A   = [ptx(id)';pty(id)'];
    naa = null(A);
    Pmat = [naa';A];
    
    Mmat = full(Smatrix'*Smatrix);
    Qmat = zeros(dim*num);
    for i=id(1:dim)
        Qmat(i,i) = 1;
    end
    Qmat(nid,nid) = -Mmat(nid,nid);
    Qmat(id,nid)  = -Pmat*Mmat(id,nid);
    Qmat(nid,id(dim+1:2*dim)) = -Mmat(nid,id)*Pmat(dim+1:2*dim,:)';
    Qmat(id,id(dim+1:2*dim))  = -Pmat*Mmat(id,id)*Pmat(dim+1:2*dim,:)';
    fext = Mmat*pert;
    fext(id) = Pmat*fext(id);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(id(1:dim)) = 0; % set force to zero
    disp(id) = Pmat'*disp(id);
    disp = disp+pert;
else
    i = 1;
    id1 = odu(i);
    while ismember(id1,in1)
        i = i+1;
        id1 = odu(i);
    end
    i = 1;
    id2 = odl(i);
    while ismember(id2,in1)
        i = i+1;
        id2 = odl(i);
    end
    id  = [2*id1-1,2*id1,2*id2-1,2*id2];
    i = 1;
    it1 = odu(i);
    while ismember(it1,in2)
        i = i+1;
        it1 = odu(i);
    end
    i = 1;
    it2 = odl(i);
    while ismember(it2,in2)
        i = i+1;
        it2 = odl(i);
    end
    it  = [2*it1-1,2*it1,2*it2-1,2*it2];
    idt = [id,it];
    nidt= setdiff(1:num*dim,idt);
    
    pert = zeros(num*dim,1);
    pert(id(1)) = -1; pert(id(2)) = 0; pert(id(3)) = 1; pert(id(4)) = 0;
    pert(it(1)) = 1; pert(it(2)) = 0; pert(it(3)) = -1; pert(it(4)) = 0;
    ptx = zeros(num*dim,1);
    pty = zeros(num*dim,1);
    ptx(id(1:2:end)) = 1/sqrt(2);
    pty(id(2:2:end)) = 1/sqrt(2);
    A   = [ptx(id)';pty(id)'];
    naa = null(A);
    Pmat = [naa';A];
    ptx = zeros(num*dim,1);
    pty = zeros(num*dim,1);
    ptx(it(1:2:end)) = 1/sqrt(2);
    pty(it(2:2:end)) = 1/sqrt(2);
    A   = [ptx(it)';pty(it)'];
    naa = null(A);
    Tmat = [naa';A];
    
    PTmat = zeros(size(Pmat)+size(Tmat));
    PTmat(1:2*dim,1:2*dim) = Pmat;
    PTmat(2*dim+(1:2*dim),2*dim+(1:2*dim)) = Tmat;
    Mmat = full(Smatrix'*Smatrix);
    Qmat = zeros(dim*num);
    for i=[id(1:dim),it(1:dim)]
        Qmat(i,i) = 1;
    end
    Qmat(nidt,nidt) = -Mmat(nidt,nidt);
    Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
    Qmat(nidt,[id(dim+1:2*dim),it(dim+1:2*dim)]) = -Mmat(nidt,idt)*PTmat([dim+1:2*dim,3*dim+1:4*dim],:)';
    Qmat(idt,[id(dim+1:2*dim),it(dim+1:2*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([dim+1:2*dim,3*dim+1:4*dim],:)';
    fext = Mmat*pert;
    fext(idt) = PTmat*fext(idt);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp([id(1:dim),it(1:dim)]) = 0;
    disp(idt) = PTmat'*disp(idt);
    disp = disp+pert;
end
fext = Mmat*pert;
ener = 0.5*disp'*Mmat*disp;
disp(1:2:end) = disp(1:2:end) - mean(disp(1:2:end));
disp(2:2:end) = disp(2:2:end) - mean(disp(2:2:end));
%{
xcen = mean(pos(:,1));
ycen = mean(pos(:,2));
rot = zeros(dim*num,1);
for i=1:num
    rot(dim*i-1) = -(pos(i,2)-ycen);
    rot(dim*i)   = pos(i,1)-xcen;
end
rot = rot/norm(rot);
disp = disp-disp'*rot*disp;
%}