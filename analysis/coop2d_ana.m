
allos = 0;
warning off;
eps   = 1e-10;

spf   = nsp/nnb;
%msgtmp = zeros(nnb,1);
%flptmp = zeros(nnb,1);
mdtemp = zeros(num,dim);
magtmp = zeros(num,1);
modtmp = zeros(num,1);
bktemp = zeros(num,1);
shtemp = zeros(num,1);
bftemp = zeros(num,1);
sbftmp = zeros(num,1);
sbtemp = zeros(num,1);
sstemp = zeros(num,1);
wttemp = zeros(num*dim,num);  % energy weight of each mode on each particle
%shistt = zeros(length(xx),1);
%    figure; hold all;
nc   = 0;

msgtmp = mean(config(81:end,:))';
for cc=82:2:100
%cc = 51;
%%
state  = config(cc,:);   % take a reference configuration
%{
pMmat  = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
comppcost;
pcost0 = pcost;   % cooperative energy of a reference configuration
tic;
for ii=1:nnb
    state1 = state;
    state1(ii) = mod(state(ii)+1,2);
    pMmat = kweak*MmatW+Smatrix'*diag(state1)*Smatrix;
    comppcost;
    flptmp(ii) = flptmp(ii)+pcost-pcost0;  % energy difference
end
toc;
%}
%%
pMmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
[V,D] = eig(full(pMmat));
[bf,stb]  = compBfactor(pMmat,[posx;posy],bondnb);  % strain b factor
% compute response field
if allos == 1  %allosteric cost
    Qmat   = zeros(dim*num);
    for i=idd
        Qmat(i,i) = 1;
    end
    Qmat(:,nid) = -pMmat(:,nid);
    fext = pMmat*pert;
    disp = linsolve(Qmat,fext);  %kweak>0
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext;
    end
    disp(idd) = pert(idd);
else % cooperative cost
    % displacement with both stimulus and target
    pQmat   = zeros(dim*num);
    %for i=itt
    %    pQmat(i,i) = 1;
    %end
    for i=idd
        pQmat(i,i) = 1;
    end
    pQmat(:,nid) = -pMmat(:,nid);  %nidt
    fext = pMmat*(pert);  %+targ
    disp = linsolve(pQmat,fext);  %pQmat\fext; %
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(pQmat)*fext; %
    end
    fx  = linsolve(pQmat,pMmat*ptx);  % tx, ty, rot
    fy  = linsolve(pQmat,pMmat*pty);
    fr  = linsolve(pQmat,pMmat*prot);
    Mm(1,1) = ptx'*fx;  Mm(1,2) = ptx'*fy;  Mm(1,3) = ptx'*fr;
    Mm(2,1) = pty'*fx;  Mm(2,2) = pty'*fy;  Mm(2,3) = pty'*fr;
    Mm(3,1) = prot'*fx; Mm(3,2) = prot'*fy; Mm(3,3) = prot'*fr;
    vb(1) = -ptx'*disp;  vb(2) = -pty'*disp;  vb(3) = -prot'*disp;
    pt    = linsolve(Mm,vb);
    
    disp = disp+pt(1)*fx+pt(2)*fy+pt(3)*fr; %...
    %+tconst(1)*tfx+tconst(2)*tfy+tconst(3)*tfr;
    disp(idd) = pert(idd)+pt(1)*ptx(idd)+pt(2)*pty(idd)+pt(3)*prot(idd);
    %disp(itt) = targ(itt)+pt(1)*tx(itt)+pt(2)*ty(itt)+pt(3)*rot(itt);
end
if nnb<3*num-2*xl
    opt.bd = 1;
else
    opt.bd = 0;
end
displ = rmvtsrot(disp,[posx;posy],opt); % displacement field with translations and rotations removed
qq  = abs(V'*displ)/norm(displ); % overlap on different modes
[qq1,ord] = sort(qq,'descend');   % rank the coupled modes
oo = num*sum((V(1:dim:end,:).^2+V(2:dim:end,:).^2).^2)';
v0 = V(:,ord(qq1>0.1));
[Vm,Dm]=eig(v0'*MmatW*v0);
%[~,id] = min(diag(Dm));
vv = v0*Vm; %(:,id);
q1 = abs(vv'*displ)/norm(displ);
[~,od1]  = sort(q1,'descend');
modtmp = modtmp+sqrt(num*sum(reshape(vv(:,od1(1)),dim,num).^2))';
tic;
for ii = 1:num*dim
    forc = diag(state)*Smatrix*V(:,ii);
    for i=1:num
        flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
        wttemp(ii,i) = wttemp(ii,i)+forc(flg)'*forc(flg)/2;
    end
end
toc;
stres = diag(diag(state)*Smatrix*displ);
sss = cell(num,1);
for i = 1:num
    smm  = Smatrix(:,(i-1)*dim+(1:dim));
    sss{i} = smm'*stres*smm; %sum(abs(stres(bid)))/2;
end
[bks,shs] = strainEnergy(sss);
[strn,vor] = computeStrain(displ,[posx;posy],bondnb);
[bk,sh] = strainEnergy(strn);
bftemp = bftemp+log10(bf);   % b-factor
sbftmp = sbftmp+log10(stb);  % strain b-factor
bktemp = bktemp+bk;   % bulk energy
shtemp = shtemp+sh;   % shear energy
sbtemp = sbtemp+bks;    % bulk stress
sstemp = sstemp+shs;   % shear stress
mdtemp = mdtemp+reshape(displ,dim,num)';
magtmp = magtmp+sqrt(sum(reshape(displ,dim,num).^2))';
%[cc,~] = histc(sh,xx);
%shistt = shistt+cc;
nc = nc+1;   % count the number of configurations averaged
end
%%
msig  = msgtmp;   % average occupation of a link
flip  = -flptmp/nc;  % average flipping cost of a link
msigc   = min(max(msig,eps),1-eps);
conserv = msigc.*log(msigc/spf)+(1-msigc).*log((1-msigc)/(1-spf));  % conservation on each link
cons  = zeros(num,1);  % conservation on each node
fcst  = zeros(num,1);  % flipping cost on each node
for i=1:num   % average link values to each node
    flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
    cons(i) = mean(conserv(flg));
    fcst(i) = mean(flip(flg));
end
wttemp = wttemp/nc;
mcst  = wttemp*fcst./sum(wttemp,2);
mdisp = mdtemp/nc;
magd  = magtmp/nc;
mode  = modtmp/nc;
bulke = bktemp/nc;
shere = shtemp/nc;
bfct  = bftemp/nc;
sbfc  = sbftmp/nc;
%shht  = shistt/nc;
strb  = sbtemp/nc;
strs  = sstemp/nc;