function [pcost,displ] = comppcost(Mmat,pert,targ,allos,pos,ids,mats)
num  = size(pos,2);
dim  = size(pos,1);
idd  = ids{1};
nid  = ids{2};
itt  = ids{3};
nit  = ids{4};
idt  = ids{5};
nidt = ids{6};
Pmat = mats{1};
Tmat = mats{2};
PTmat = mats{3};
%pMmat = kweak*MmatW+Smat'*diag(nstat)*Smat;
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
    diff = disp(itt)-targ(itt);
    pcost = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
    displ = rmvtsrot(disp,pos);
else
    % displacement with stimulus
    Qmat = zeros(dim*num);
    for i=idd(1:2*dim+1)
        Qmat(i,i) = 1;
    end
    Qmat(nid,nid) = -Mmat(nid,nid);
    Qmat(idd,nid) = -Pmat*Mmat(idd,nid);
    Qmat(nid,idd(3*dim:4*dim)) = -Mmat(nid,idd)*Pmat(3*dim:4*dim,:)';
    Qmat(idd,idd(3*dim:4*dim)) = -Pmat*Mmat(idd,idd)*Pmat(3*dim:4*dim,:)';
    fext = Mmat*pert;
    fext(idd) = Pmat*fext(idd);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    ff = disp(idd(1:2*dim+1)); 
    disp(idd(1:2*dim+1)) = 0; % set force to zero
    disp(idd) = Pmat'*disp(idd);
    disp = disp+pert;
    displ = rmvtsrot(disp,pos);
    pes  = 0.5*disp'*Mmat*disp;
    
    % energy of only targ
    Qmat = zeros(dim*num);
    for i=itt(1:2*dim+1)
        Qmat(i,i) = 1;
    end
    Qmat(nit,nit) = -Mmat(nit,nit);
    Qmat(itt,nit) = -Tmat*Mmat(itt,nit);
    Qmat(nit,itt(3*dim:4*dim)) = -Mmat(nit,itt)*Tmat(3*dim:4*dim,:)';
    Qmat(itt,itt(3*dim:4*dim)) = -Tmat*Mmat(itt,itt)*Tmat(3*dim:4*dim,:)';
    fext = Mmat*targ;
    fext(itt) = Tmat*fext(itt);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp(itt(1:2*dim+1)) = 0;
    disp(itt) = Tmat'*disp(itt);
    disp = disp+targ;
    pet  = 0.5*disp'*Mmat*disp;
    
    % displacement with both stimulus and target
    Qmat = zeros(dim*num);
    for i=[idd(1:2*dim+1),itt(1:2*dim+1)]
        Qmat(i,i) = 1;
    end
    Qmat(nidt,nidt) = -Mmat(nidt,nidt);
    Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
    Qmat(nidt,[idd(3*dim:4*dim),itt(3*dim:4*dim)]) = -Mmat(nidt,idt)*PTmat([3*dim:4*dim,7*dim:8*dim],:)';
    Qmat(idt,[idd(3*dim:4*dim),itt(3*dim:4*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([3*dim:4*dim,7*dim:8*dim],:)';
    fext = Mmat*(pert+targ);
    fext(idt) = PTmat*fext(idt);
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext; %
    end
    disp([idd(1:2*dim+1),itt(1:2*dim+1)]) = 0;
    disp(idt) = PTmat'*disp(idt);
    disp = disp+pert+targ;
    pest = 0.5*disp'*Mmat*disp;
    pcost = pes+pet-pest; %cost
end