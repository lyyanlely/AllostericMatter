%% read the numerical data
dirc  = '../Allostery/Nov27/';  %Nov18
ename = 'Costs';
zname = 'MeanZ';
pname = 'PartRatio';
cname = 'Config';
plnam = 'PrLocal';
dpnam = 'DispMap';
zlnam = 'Zlocal';
rsnam = 'Response';
bfnam = 'Bfactor';
sfnam = 'Sfactor';
ssnam = 'ShearSpring';
sbnam = 'ShearBond';
bsnam = 'BulkSpring';
bbnam = 'BulkBond';
ncnam = 'ChangeCost';
esnam = 'Stimulus';
etnam = 'Target';
stnam = 'Cooperativity';
cbnam = 'CorrelBond';
csnam = 'CorrelSpring';
dtype = '.dat';
dd    = 1;
dt    = 1;
kweak = 0.0001;
dname = sprintf('%.2f',dd);
tname = sprintf('%.2f',dt);
kname = sprintf('%.4f',kweak);

%npara = 1;
nsp = 1000;
rp    = 50;
%{
beta_ini=1000;
beta  = zeros(1,npara);
if npara>3
    b0   = beta_ini^(1/(npara-2));
    if b0>1
        beta(2:npara-1) = b0.^(0:npara-3);
    else
        beta(2:npara-1) = b0.^(3-npara:0);
    end
    beta(npara) = 1e12;
elseif npara==3
    beta = [1,beta_ini,1e12];
elseif npara==2
    beta = [beta_ini,1e12];
else
    beta = 1e12;
end
%}
%beta = [0.04,0.05,0.06,0.065,0.07,0.075,0.08,0.085,0.088,0.089,0.09,0.091,0.092,0.095,0.1,0.11,0.12,0.14,0.2,0.4];
beta = 0.005; %[0.001:0.001:0.009,0.01:0.01:0.1,0.11,0.12,0.14,0.16,0.2,0.25,0.3]; % [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3]; %[0.05,0.4]; %[0,0.04:0.015:0.16]; %

nx  = 20; %[4,6,8,10,12];
ny  = 20; %[4,6,8,10,12]; %[4,6,8,10,12];%[3,4,5,6,7,8,9,10,11,12];

lx  = length(nx);
npara = length(beta);

pos = [posx;posy];
ids = {idd,nid,itt,nit,idt,nidt};
mats= {Pmat,Tmat,PTmat};
%% compute the thermodynamic factors
allos = 1;
warning off;
eps   = 1e-10;

xname  = sprintf('%d',nx);
yname  = sprintf('%d',ny);
num    = ny*nx;
%nsp    = round(num*5/2); %(2*nx*ny(j)-4):(3*nx*ny(j)-2*nx-1);
    
spf    = nsp/nnb;
nname  = sprintf('%d',nsp);

msgtmp = zeros(nnb,1);
flptmp = zeros(nnb,1);
mdtemp = zeros(num,dim);
magtmp = zeros(num,1);
modtmp = zeros(num,1);
bktemp = zeros(num,1);
shtemp = zeros(num,1);
bftemp = zeros(num,1);
sbftmp = zeros(num,1);
sbtemp = zeros(num,1);
sstemp = zeros(num,1);
zltmp  = zeros(1,num);
wttemp = zeros(num*dim,num);  % energy weight of each mode on each particle
%shistt = zeros(length(xx),1);
%    figure; hold all;
nc   = 0;
tnam = sprintf('%.4f',beta);

for r = 1:rp
    rname = sprintf('%03d',r);
    filename = [dirc,ename,'_',xname,'_',yname,'_',nname,'_',tnam,'_',kname,'_',rname,dtype];  %,'_',tnam
    zlocname = [dirc,zlnam,'_',xname,'_',yname,'_',nname,'_',tnam,'_',kname,'_',rname,dtype];
    dispname = [dirc,dpnam,'_',xname,'_',yname,'_',nname,'_',tnam,'_',kname,'_',rname,dtype];
    confname = [dirc,cname,'_',xname,'_',yname,'_',nname,'_1.00_1.00_',tnam,'_',kname,'_1_',rname,dtype]; %
    if exist(confname,'file')
        s = dir(confname);
        if s.bytes>1
            %{
            cdata = dlmread(filename);
            if min(size(cdata))<npara
                cdata = reshape(cdata,[],1);
            else
                cdata = cdata';
            end
            zldat = dlmread(zlocname);
            if min(size(zldat))<npara
                zldat = reshape(zldat,1,[]);
            end
            % mean displacement
            dpdat = dlmread(dispname);
            if min(size(dpdat))<npara
                dpdat = reshape(dpdat,1,[]);
            end
            %}
            config1 = [dlmread(confname)';ones(nb1,100)]';
            if nc==0
                config = config1(51:end,:);
                nc = 1;
            else
                config = [config;config1(51:end,:)];
            end
        end
    end
    %
    msgtmp = msgtmp+mean(config(51:end,:))';
    for cc=60:10:100
        %%
        state  = config(cc,:);   % take a reference configuration
        % cost of flipping a bond (comment out this part if you don't need,
        % very slow)
        pMmat  = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
        [pcost,displ] = comppcost(Mmat,pert,targ,allos,pos,ids,mats);
        pcost0 = pcost;   % cooperative energy of a reference configuration
        tic;
        for ii=1:nnb
            state1 = state;
            state1(ii) = mod(state(ii)+1,2);
            pMmat = kweak*MmatW+Smatrix'*diag(state1)*Smatrix;
            [pcost,~] = comppcost(Mmat,pert,targ,allos,pos,ids,mats);
            flptmp(ii) = flptmp(ii)+pcost-pcost0;  % energy difference
        end
        toc;
        %
        %%
        pMmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
        [V,D] = eig(full(pMmat));
        % b-factor and strain B-factor
        [bf,stb]  = compBfactor(pMmat,[mod(posx+8,nx);posy],bondnb);  % strain b factor
        %[pcost,displ] = comppcost(Mmat,pert,targ,allos,pos,ids,mats);
        qq  = abs(V'*displ)/norm(displ); % overlap on different modes
        [qq1,ord] = sort(qq,'descend');   % rank the coupled modes
        oo = num*sum(V.^4)';
        v0 = V(:,ord(qq1>0.1));
        [Vm,Dm]=eig(v0'*MmatW*v0);
        %[~,id] = min(diag(Dm));
        vv = v0*Vm; %(:,id);
        q1 = abs(vv'*displ)/norm(displ);
        [~,od1]  = sort(q1,'descend');
        modtmp = modtmp+sqrt(num*sum(reshape(vv(:,od1(1)),dim,num).^2))';
        %{
        tic;
        for ii = 1:num*dim
            forc = diag(state)*Smatrix*V(:,ii);
            for i=1:num
                flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
                wttemp(ii,i) = wttemp(ii,i)+forc(flg)'*forc(flg)/2;
            end
        end
        toc;
        %}
        % compute pseudo-strain
        stres = diag(diag(state)*Smatrix*displ);
        sss = cell(num,1);
        for i = 1:num
            smm  = Smatrix(:,(i-1)*dim+(1:dim));
            sss{i} = smm'*stres*smm; %sum(abs(stres(bid)))/2;
        end
        [strn,vor] = computeStrain(displ,[mod(posx+8,nx);posy],bondnb);
        % pseudo-energy
        [bks,shs] = strainEnergy(sss);
        [bk,sh] = strainEnergy(strn);
        bftemp = bftemp+bf;   % b-factor
        sbftmp = sbftmp+stb;  % strain b-factor
        bktemp = bktemp+bk;   % bulk energy
        shtemp = shtemp+sh;   % shear energy
        sbtemp = sbtemp+bks;    % bulk stress
        sstemp = sstemp+shs;   % shear stress
        mdtemp = mdtemp+reshape(displ,dim,num)';
        magtmp = magtmp+sqrt(sum(reshape(displ,dim,num).^2))';
        zltmp  = zltmp+zldat;  %(n,:)
        %[cc,~] = histc(sh,xx);
        %shistt = shistt+cc;
        nc = nc+1;   % count the number of configurations averaged
    end
    %}
end
%%
%
msig  = msgtmp/rp;   % average occupation of a link
flip  = -flptmp;  % average flipping cost of a link
msigc   = min(max(msig,eps),1-eps);
conserv = msigc.*log(msigc/spf)+(1-msigc).*log((1-msigc)/(1-spf));  % conservation on each link
cons  = zeros(num,1);  % conservation on each node
% comment out if not computing fcst
fcst  = zeros(num,1);  % flipping cost on each node
for i=1:num   % average link values to each node
    flg = ismember(bondnb(:,1),i)|ismember(bondnb(:,2),i);
    cons(i) = mean(conserv(flg));
    fcst(i) = mean(flip(flg));
end
%wttemp = wttemp/nc;
%mcst  = wttemp*fcst./sum(wttemp,2);
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

zloc  = zltmp/nc;
%}