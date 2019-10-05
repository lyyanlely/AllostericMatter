%% read the numerical data
dirc  = '../Allostery/3D/';
cname = 'Costs3D';
sname = 'Config3D';
cdnam = 'Zlocal3D';
esnam = 'Stimulus3D';
etnam = 'Target3D';
stnam = 'Cooperativity3D';
allos = 1;
spf   = 0.8;
kweak = 0.0001;
eps   = 1e-10;
kname = sprintf('%.4f',kweak);
aname = sprintf('%d',allos);
warning off

rp    = 20;
beta = 0.004; %[0.001:0.001:0.009,0.01:0.01:0.1,0.11,0.12,0.14,0.16,0.2,0.25,0.3]; % [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3];

nx  = 6; %[4,6,8,10,12];
ny  = 12; %[4,6,8,10,12]; 
nz  = 12;

xx  = 10.^(-4:0.5:2);

lx  = length(nx);
npara = length(beta);
%% compute the thermodynamic factors
num  = zeros(lx,1);
nsp  = zeros(lx,1);
cost = cell(lx,1);
hcap = cell(lx,1);
estim= cell(lx,1);   % energy from stimulus
etarg= cell(lx,1);   % energy from target
estig= cell(lx,1);   % energy from stimulus+target

msig = cell(lx,npara); % mean occupancy
flip = cell(lx,npara); % flip cost
cons = cell(lx,npara); % conservation
fcst = cell(lx,npara); % cooperative energy cost to flip state
mcst = cell(lx,npara); % cooperative energy cost on each mode
zloc = cell(lx,npara); % coordination distribution
mdisp= cell(lx,npara); % mean displacement
magd = cell(lx,npara); % magnitude of displacement
mode = cell(lx,npara); % magnitude of the best fit mode
bulke= cell(lx,npara); % bulk energy
shere= cell(lx,npara); % shear energy
bfct = cell(lx,npara); % b-factor
sbfc = cell(lx,npara); % strain b-factor
strb = cell(lx,npara); % stress
strs = cell(lx,npara); % stress
shht = cell(lx,npara); % histgram
ovlp = cell(rp,10);   % overlap

qaa = 0;
qab = 0;

pp = 0;
for j = 1:lx
    xname  = sprintf('%d',nx(j));
    yname  = sprintf('%d',ny(j));
    zname  = sprintf('%d',nz(j));
    num(j) = nz(j)*ny(j)*nx(j);
    nlnk   = nnb;
    nsp(j) = round(nlnk*spf);
    ln  = 1;
    cost{j} = zeros(ln,npara);
    estim{j}= zeros(ln,npara);
    etarg{j}= zeros(ln,npara);
    estig{j}= zeros(ln,npara);
    zloc{j} = zeros(npara,num(j));
    
    nname  = sprintf('%d',nsp(j)); %/(3*num(j)-2*nx(j))
    ctemp  = zeros(1,npara);
    htemp  = zeros(1,npara);
    estmp  = zeros(1,npara);
    ettmp  = zeros(1,npara);
    sttmp  = zeros(1,npara);
    zltmp  = zeros(npara,num(j));
    nr     = zeros(npara,1);
    for n = 1:npara
        tnam = sprintf('%.4f',beta(n));
        msgtmp = zeros(nnb,1);
        flptmp = zeros(nnb,1);
        mdtemp = zeros(num(j),dim);
        magtmp = zeros(num(j),1);
        modtmp = zeros(num(j),1);
        bktemp = zeros(num(j),1);
        shtemp = zeros(num(j),1);
        bftemp = zeros(num(j),1);
        sbftmp = zeros(num(j),1);
        sbtemp = zeros(num(j),1);
        sstemp = zeros(num(j),1);
        wttemp = zeros(num(j)*dim,num(j));  % energy weight of each mode on each particle
        shistt = zeros(length(xx),1);
        %    figure; hold all;
        nc   = 0;
        %%
        for r = 1:rp
            rname = sprintf('%03d',r);
            costname = [dirc,cname '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
            mzname = [dirc,cdnam '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
            confname = [dirc,sname '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
            if allos~=1
                esname = [dirc,esnam '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
                etname = [dirc,etnam '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
                stname = [dirc,stnam '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
            end
            %{
            if exist(costname,'file')
                s = dir(costname);
                if s.bytes>1
                    cdata = dlmread(costname);
                    %plot(1000:1000:100000,cdata);
                    zldat = dlmread(mzname);
                    if nr(n)==0
                        ztemp = sum(zldat(71:100,:));
                    else
                        ztemp = ztemp+sum(zldat(71:100,:));
                    end
                    % energy on stimulus target and both
                    %
                    if allos~=1
                        esdat = dlmread(esname);
                        etdat = dlmread(etname);
                        stdat = dlmread(stname);
                    end
                    %
                    
                    cl    = size(cdata);
                    cl    = cl(2);
                    ctemp(n) = ctemp(n)+mean(cdata(round(cl/2)+1:cl));  %(n)
                    htemp(n) = htemp(n)+var(cdata(round(cl/2)+1:cl));
                    %zltmp(n,:) = zltmp(n,:)+zldat;  %(n,:)
                    %
                    if allos~=1
                        estmp(n) = estmp(n)+mean(esdat(round(cl/2)+1:cl));
                        ettmp(n) = ettmp(n)+mean(etdat(round(cl/2)+1:cl));
                        sttmp(n) = sttmp(n)+mean(stdat(round(cl/2)+1:cl));
                    end
                    %
                    nr(n)       = nr(n)+1;
                end
            end
            %}
            if exist(confname,'file')
                config1 = dlmread(confname);
                %qaa = qaa + (config(17,:)*config(50,:)'/spf/nlnk-spf)/(1-spf);
                
                %if nc>1
                %    qab = qab + (config(50,:)*config1(50,:)'/spf/nlnk-spf)/(1-spf);
                %end
                if nc==0
                    config = config1;
                    nc = 1;
                else
                    config = [config;config1];
                end
                %msgtmp = msgtmp+mean(config(21:50,:))';
                %% compute the single flipping cost
                %{
                tic;
                %i = 50;
                for cc = 23:3:50
                state  = config(cc,:);
                %
                Mmat  = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
                compcost3D;
                pcost0 = pcost;
                displ = rmvtsrot(disp,pos);
                pp = num(j)*sum((displ(1:dim:end).^2+displ(2:dim:end).^2+displ(3:dim:end).^2).^2)/sum(displ.^2)^2;
                %vMat = xlaplacian(Mmat,15);
                %q0 = abs(vMat'*displ)/norm(displ);
                %[~,id] = max(q0);
                [vMat,~] = eigs(Mmat,10+6,'sm');
                q0 = abs(vMat'*displ)/norm(displ);
                %qq = num(j)*sum((vMat(1:dim:end,id).^2+vMat(2:dim:end,id).^2+vMat(3:dim:end,id).^2).^2);
                ovlp{r,(cc-20)/3} = q0.^2;
                %
                tic;
                for ii=1:nnb
                    state1 = state;
                    state1(ii) = mod(state(ii)+1,2);
                    Mmat = kweak*MmatW+Smatrix'*diag(state1)*Smatrix;
                    compcost3D;
                    flptmp(ii) = flptmp(ii)+pcost-pcost0;
                end
                toc;
                %
                %%
                %{
                Mmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
                [V,D] = eig(full(Mmat));
                [bf,stb]  = compBfactor(Mmat,pos,bondnb);  % strain b factor
                % compute response field
                if allos == 1  %allosteric cost
                    Qmat   = zeros(dim*num(j));
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
                else % cooperative cost
                    % displacement with both stimulus and target
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
                end
                displ = rmvtsrot(disp,pos); % displacement field with translations and rotations removed
                qq  = abs(V'*displ)/norm(displ); % overlap on different modes
                [qq1,ord] = sort(qq,'descend');   % rank the coupled modes
                oo = num*sum((V(1:3:end,:).^2+V(2:3:end,:).^2+V(3:3:end,:).^2).^2)';
                v0 = V(:,ord(qq1>0.1));
                [Vm,Dm]=eig(v0'*MmatW*v0);
                %[~,id] = min(diag(Dm));
                vv = v0*Vm; %(:,id);
                q1 = abs(vv'*displ)/norm(displ);
                [~,od1]  = sort(q1,'descend');
                modtmp = modtmp+sqrt(num*sum(reshape(vv(:,od1(1)),dim,num).^2))';
                %{
                tic;
                for ii = 1:num(j)*dim
                    forc = diag(state)*Smatrix*V(:,ii);
                    for i=1:num(j)
                        flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
                        wttemp(ii,i) = wttemp(ii,i)+forc(flg)'*forc(flg)/2;
                    end
                end
                toc;
                %}
                %{
                stres = diag(state)*Smatrix*displ;
                ss = zeros(num(j),1);
                for i = 1:num(j)
                    bid = bondnb(:,1)==i|bondnb(:,2)==i;
                    ss(i) = sum(abs(stres(bid)))/2;
                end
                %}
                stres = diag(diag(state)*Smatrix*displ);
                sss = cell(num(j),1);
                for i = 1:num(j)
                    %bid = bondnb(:,1)==i|bondnb(:,2)==i;
                    smm  = Smatrix(:,(i-1)*dim+(1:dim));
                    sss{i} = smm'*stres*smm; %sum(abs(stres(bid)))/2;
                end
                [bks,shs] = strainEnergy(sss);
                [strn,vor] = computeStrain(displ,pos,bondnb);
                [bk,sh] = strainEnergy(strn);
                bftemp = bftemp+log10(bf);   % b-factor
                sbftmp = sbftmp+log10(stb);  % strain b-factor
                bktemp = bktemp+bk;   % bulk energy
                shtemp = shtemp+sh;   % shear energy
                sbtemp = sbtemp+bks;    % bulk stress
                sstemp = sstemp+shs;   % shear stress
                %vortmp = vortmp+vor;  % vorticity
                mdtemp = mdtemp+reshape(displ,dim,num(j))';
                magtmp = magtmp+sqrt(sum(reshape(displ,dim,num).^2))';
                %[cc,~] = histc(sh,xx);
                %shistt = shistt+cc;
            %}
                nc = nc+1;
                end
                %}
                %toc;
            end
        end
        %%
        %{
        msig{j,n}  = msgtmp;
        flip{j,n}  = -flptmp/nc;
        msigc   = min(max(msig{j,n},eps),1-eps);
        conserv = msigc.*log(msigc/spf)+(1-msigc).*log((1-msigc)/(1-spf));
        cons{j,n}  = zeros(num(j),1);
        fcst{j,n}  = zeros(num(j),1);
        for i=1:num(j)
            flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
            cons{j,n}(i) = mean(conserv(flg));
            fcst{j,n}(i) = mean(flip{j,n}(flg));
        end
        wttemp = wttemp/nc;
        mcst{j,n}  = wttemp*fcst{j,n}./sum(wttemp,2);
        mdisp{j,n} = mdtemp/nc;
        magd{j,n}  = magtmp/nc;
        mode{j,n}  = modtmp/nc;
        bulke{j,n} = bktemp/nc;
        shere{j,n} = shtemp/nc;
        bfct{j,n}  = bftemp/nc;
        sbfc{j,n}  = sbftmp/nc;
        shht{j,n}  = shistt/nc;
        strb{j,n}  = sbtemp/nc;
        strs{j,n}  = sstemp/nc;
        %}
        %{
        ztemp = zeros(num,1);
        for n0 = 1:nnb
            ztemp(bondnb(n0,1)) = ztemp(bondnb(n0,1))+msig{j,n}(n0);
            ztemp(bondnb(n0,2)) = ztemp(bondnb(n0,2))+msig{j,n}(n0);
        end
        %}
        %zloc{j,n}  = ztemp/nr(n)/30;
        %%
    end
    cost{j}(1,:) = ctemp./nr';
    hcap{j}(1,:) = htemp./nr';
    if allos~=1
        estim{j}(1,:) = estmp./nr';
        etarg{j}(1,:) = ettmp./nr';
        estig{j}(1,:) = sttmp./nr';
    end
    %zloc{j} = zltmp./repmat(nr,1,num(j));
    
end