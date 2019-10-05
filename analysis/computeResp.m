function displ = computeResp(Mmat,pos,allos,pert,idd,Pmat)
num = size(pos,2);
dim = size(pos,1);
nid = setxor(1:num*dim,idd);
%{
idx = find(pert~=0);
id  = union(ceil(idx/dim),[])';
idd = zeros(1,dim*length(id));
for i = 1:length(id)
    idd((i-1)*dim+1:i*dim) = (id(i)-1)*dim+1:id(i)*dim;
end
nid = setxor(1:dim*num,idd);
%}

if allos==1
    Qmat   = zeros(dim*num);
    for i=idd
        Qmat(i,i) = 1;
    end
    Qmat(:,nid) = -Mmat(:,nid);
    fext = Mmat*pert;
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext;
    end
    disp(idd) = pert(idd);
    displ = rmvtsrot(disp,pos);
else
    Qmat = zeros(dim*num);
    for i=idd(1:2*dim)
        Qmat(i,i) = 1;
    end
    Qmat(nid,nid) = -Mmat(nid,nid);
    Qmat(idd,nid) = -Pmat*Mmat(idd,nid);
    Qmat(nid,idd(2*dim+1:4*dim)) = -Mmat(nid,idd)*Pmat(2*dim+1:4*dim,:)';
    Qmat(idd,idd(2*dim+1:4*dim)) = -Pmat*Mmat(idd,idd)*Pmat(2*dim+1:4*dim,:)';
    fext = Mmat*pert;
    fext(idd) = Pmat*fext(idd);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(idd(1:2*dim)) = 0; % set force to zero
    disp(idd) = Pmat'*disp(idd);
    disp = disp+pert;
    displ = rmvtsrot(disp,pos);
end