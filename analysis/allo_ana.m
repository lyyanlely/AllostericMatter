%% read the numerical data
dirc  = '../Allostery/Nov27/';
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
rp    = 20;
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
beta = 0.01; %[0.001:0.001:0.009,0.01:0.01:0.1,0.11,0.12,0.14,0.16,0.2,0.25,0.3]; % [0.0001,0.0003,0.001,0.003,0.01,0.03,0.1,0.3]; %[0.05,0.4]; %[0,0.04:0.015:0.16]; %

nx  = 20; %[4,6,8,10,12];
ny  = 20; %[4,6,8,10,12]; %[4,6,8,10,12];%[3,4,5,6,7,8,9,10,11,12];

lx  = length(nx);
npara = length(beta);
%% compute the thermodynamic factors
num  = zeros(lx,1);
nsp  = cell(lx,1);
cost = cell(lx,1);
hcap = cell(lx,1);
cdz  = cell(lx,1);
part = cell(lx,1);
estim= cell(lx,1);
etarg= cell(lx,1);
estig= cell(lx,1);
zloc = cell(lx,1);
ploc = cell(lx,1);
dpmp = cell(lx,1);
resp = cell(lx,1);
bfac = cell(lx,1);
sfac = cell(lx,1);
shsp = cell(lx,1);
shbd = cell(lx,1);
bksp = cell(lx,1);
bkbd = cell(lx,1);
crr  = cell(lx,1);  % distance r
crbd = cell(lx,1);  % cumulative probability over distance for bond
crsp = cell(lx,1);  % ... for spring
dcb  = cell(lx,1);
dcs  = cell(lx,1);  % dimension
flct = cell(lx,1);
ncst = cell(lx,1);

for j = 1:lx
    xname  = sprintf('%d',nx(j));
    yname  = sprintf('%d',ny(j));
    num(j) = ny(j)*nx(j);
    nsp{j} = round(num(j)*5/2); %(2*nx*ny(j)-4):(3*nx*ny(j)-2*nx-1);
    nlnk   = 1102; %3*num(j)-2*nx(j);
    ln  = length(nsp{j});
    cost{j} = zeros(ln,npara);
    cdz{j}  = zeros(ln,npara);
    part{j} = zeros(ln,npara);
    estim{j}= zeros(ln,npara);
    etarg{j}= zeros(ln,npara);
    estig{j}= zeros(ln,npara);
    zloc{j} = zeros(npara,num(j));
    ploc{j} = zeros(npara,num(j));
    dpmp{j} = zeros(npara,2*num(j));
    resp{j} = zeros(nx,num(j));
    bfac{j} = zeros(npara,num(j));
    sfac{j} = zeros(npara,nlnk);
    shsp{j} = zeros(npara,nlnk);
    shbd{j} = zeros(npara,nlnk);
    bksp{j} = zeros(npara,nlnk);
    bkbd{j} = zeros(npara,nlnk);
    flct{j} = zeros(npara,nlnk);
    ncst{j} = zeros(npara,num(j));
    for k = 1:ln
        nname  = sprintf('%d',nsp{j}(k)); %/(3*num(j)-2*nx(j))
        ctemp  = zeros(1,npara);
        htemp  = zeros(1,npara);
        ztemp  = zeros(1,npara);
        ptemp  = zeros(1,npara);
        estmp  = zeros(1,npara);
        ettmp  = zeros(1,npara);
        sttmp  = zeros(1,npara);
        zltmp  = zeros(npara,num(j));
        pltmp  = zeros(npara,num(j));
        dptmp  = zeros(npara,2*num(j));
        rstmp  = zeros(nx,num(j));
        bftmp  = zeros(npara,num(j));
        sftmp  = zeros(npara,nlnk);
        sstmp  = zeros(npara,nlnk);
        sbtmp  = zeros(npara,nlnk);
        bstmp  = zeros(npara,nlnk);
        bbtmp  = zeros(npara,nlnk);
        nctmp  = zeros(npara,nlnk);
        nr     = zeros(npara,1);
        for n = 1:npara
            tnam = sprintf('%.4f',beta(n));
        for r = 1:rp
            rname = sprintf('%03d',r);
            filename = [dirc,ename,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];  %,'_',tnam
            cordname = [dirc,zname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            partname = [dirc,pname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            zlocname = [dirc,zlnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            prlcname = [dirc,plnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            dispname = [dirc,dpnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            respname = [dirc,rsnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            bfacname = [dirc,bfnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];            
            sfacname = [dirc,sfnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            chconame = [dirc,ncnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            estiname = [dirc,esnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            etarname = [dirc,etnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            estgname = [dirc,stnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_',rname,dtype];
            shspname = [dirc,ssnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            shbdname = [dirc,sbnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            bkspname = [dirc,bsnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            bkbdname = [dirc,bbnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            crbdname = [dirc,cbnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            crspname = [dirc,csnam '_' xname '_' yname '_' nname '_' dname '_' tname,'_',tnam,'_',kname,'_',rname,dtype];
            if exist(filename,'file')
                s = dir(filename);
                if s.bytes>1
                    cdata = dlmread(filename);
                    if min(size(cdata))<npara
                        cdata = reshape(cdata,[],1);
                    else
                        cdata = cdata';
                    end
                    zdata = dlmread(cordname);
                    if min(size(zdata))<npara
                        zdata = reshape(zdata,[],1);
                    else
                        zdata = zdata';
                    end
                    pdata = dlmread(partname);
                    if min(size(pdata))<npara
                        pdata = reshape(pdata,[],1);
                    else
                        pdata = pdata';
                    end
                    zldat = dlmread(zlocname);
                    if min(size(zldat))<npara
                        zldat = reshape(zldat,1,[]);
                    end
                    pldat = dlmread(prlcname);
                    if min(size(pldat))<npara
                        pldat = reshape(pldat,1,[]);
                    end
                    % mean displacement
                    dpdat = dlmread(dispname);
                    if min(size(dpdat))<npara
                        dpdat = reshape(dpdat,1,[]);
                    end
                    % one bond change cost
                    %ncdat = dlmread(chconame);
                    %    ncdat = reshape(ncdat,1,[]);
                    %rsdat = dlmread(respname);
                    %if min(size(rsdat))<npara
                    %    rsdat = reshape(rsdat,nx,[]);
                    %end
                    %
                    bfdat = dlmread(bfacname);
                        bfdat = reshape(bfdat,1,[]);
                    sfdat = dlmread(sfacname);
                        sfdat = reshape(sfdat,1,[]);
                    ssdat = dlmread(shspname);
                        ssdat = reshape(ssdat,1,[]);
                    sbdat = dlmread(shbdname);
                        sbdat = reshape(sbdat,1,[]);
                    bsdat = dlmread(bkspname);
                        bsdat = reshape(bsdat,1,[]);
                    bbdat = dlmread(bkbdname);
                        bbdat = reshape(bbdat,1,[]);
                    %
                    % energy on stimulus target and both
                    esdat = dlmread(estiname);
                    if min(size(esdat))<npara
                        esdat = reshape(esdat,[],1);
                    else
                        esdat = esdat';
                    end
                    etdat = dlmread(etarname);
                    if min(size(etdat))<npara
                        etdat = reshape(etdat,[],1);
                    else
                        etdat = etdat';
                    end
                    stdat = dlmread(estgname);
                    if min(size(stdat))<npara
                        stdat = reshape(stdat,[],1);
                    else
                        stdat = stdat';
                    end
                    %
                    % correlation data for fractal dimension
                    cbdat = dlmread(crbdname);
                    csdat = dlmread(crspname);
                    if nr(n) == 0 
                        crtmp = cbdat(1,:);
                        cbtmp = cbdat(2:end,:);
                        cstmp = csdat(2:end,:);
                    else
                        cbtmp = cbtmp+cbdat(2:end,:);
                        cstmp = cstmp+csdat(2:end,:);
                    end
                    %}
                    
                    cl    = size(cdata);
                    cl    = cl(1);
                    %for n = 1:npara
                        ctemp(n) = ctemp(n)+mean(cdata(round(cl/2)+1:cl,:));  %(n)
                        htemp(n) = htemp(n)+var(cdata(round(cl/2)+1:cl,:));
                        ztemp(n) = ztemp(n)+mean(zdata(round(cl/2)+1:cl,:));
                        ptemp(n) = ptemp(n)+mean(pdata(round(cl/2)+1:cl,:));
                        zltmp(n,:) = zltmp(n,:)+zldat;  %(n,:)
                        pltmp(n,:) = pltmp(n,:)+pldat;
                        dptmp(n,:) = dptmp(n,:)+dpdat;
                        %nctmp(n,:) = nctmp(n,:)+ncdat(1,:);
                        %
                        %rstmp      = rstmp+rsdat;
                        bftmp(n,:) = bftmp(n,:)+bfdat(1,:);
                        sftmp(n,:) = sftmp(n,:)+sfdat(1,:);
                        sstmp(n,:) = sstmp(n,:)+ssdat(1,:);
                        sbtmp(n,:) = sbtmp(n,:)+sbdat(1,:);
                        bstmp(n,:) = bstmp(n,:)+bsdat(1,:);
                        bbtmp(n,:) = bbtmp(n,:)+bbdat(1,:);
                        estmp(n) = estmp(n)+mean(esdat(round(cl/2)+1:cl,:));
                        ettmp(n) = ettmp(n)+mean(etdat(round(cl/2)+1:cl,:));
                        sttmp(n) = sttmp(n)+mean(stdat(round(cl/2)+1:cl,:));
                        %}
                        %if n==1
                            nr(n)       = nr(n)+1;
                        %end
                    %end
                end
            end
        end
        end
        cost{j}(k,:) = ctemp./nr';
        hcap{j}(k,:) = htemp./nr';
        cdz{j}(k,:)  = ztemp./nr';
        part{j}(k,:) = ptemp./nr';
        estim{j}(k,:) = estmp./nr';
        etarg{j}(k,:) = ettmp./nr';
        estig{j}(k,:) = sttmp./nr';
        zloc{j} = zltmp./repmat(nr,1,num(j));
        ploc{j} = pltmp./repmat(nr,1,num(j));
        dpmp{j} = dptmp./repmat(nr,1,2*num(j));
        resp{j} = rstmp./nr;
        %flct{j} = nctmp./repmat(nr,1,3*num(j)-2*nx(j))-cost{j};
        %
        crr{j}  = crtmp;
        crbd{j} = cbtmp./nr(1);
        crsp{j} = cstmp./nr(1);
        dcb{j}  = zeros(1,size(crbd{j},1));
        dcs{j}  = zeros(1,size(crsp{j},1));
        lcr = log10(crr{j});
        r0  = (lcr(1))/2;
        for ii = 1:size(crbd{j},1)
            p = polyfit(lcr,log10(crbd{j}(ii,:)),4);
            dcb{j}(ii) = 4*p(1)*r0^3+3*p(2)*r0^2+2*p(3)*r0+p(4);
            p = polyfit(lcr,log10(crsp{j}(ii,:)),4);
            dcs{j}(ii) = 4*p(1)*r0^3+3*p(2)*r0^2+2*p(3)*r0+p(4);
        end
        bfac{j} = bftmp./repmat(nr,1,num(j));
        sfac{j} = sftmp./repmat(nr,1,nlnk);
        shsp{j} = sstmp./repmat(nr,1,nlnk);
        shbd{j} = sbtmp./repmat(nr,1,nlnk);
        bksp{j} = bstmp./repmat(nr,1,nlnk);
        bkbd{j} = bbtmp./repmat(nr,1,nlnk);
        %{
        for i = 1:num(j)
            ncst{j}(:,i) = mean(flct{j}(:,nbsite{i}),2);
        end
        %}
    end
    
    % compute the eigenmodes
    %{
    stat = config(end,:,1);
    Qmat   = zeros(dim*num);
    for i=idd
        Qmat(i,i) = 1;
    end
    Mmat = kweak*MmatW+Smatrix'*diag(stat)*Smatrix;
    Qmat(:,nid) = -Mmat(:,nid);
    fext = Mmat*pert;
    disp = pinv(Qmat)*fext; %linsolve(Qmat,fext); % Qmat\fext; %
    force = zeros(1,dim*num);
    force(idd) = disp(idd);
    quiver(posx,posy,force(1:2:end),force(2:2:end),'linewidth',2,'AutoScaleFactor',1,'color',[0.8,0.7,0.])
    displ = zeros(1,dim*num);
    xid = mod(nid,2)==1;
    yid = mod(nid,2)==0;
    dx = mean(disp(nid(xid)));
    dy = mean(disp(nid(yid)));
    displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
    [V,d] = eig(Mmat);
    d=diag(d);
    omeg = sqrt(d);
    overlap=zeros(1,num*dim);
    for i=1:num*dim
        overlap(i)=displ*V(:,i);
    end
    partr = zeros(1,num*dim);
    for i=1:num*dim
        partr(i)  = 1./num/dim/sum(V(:,i).^4);
    end
    %}
end

%% map of links from Carolina
bdnb0 = floor((bondnb-1)/nx)+mod(bondnb-1,nx)*ny+1;
bdnb1 = sort(bdnb0,2);
[~,ord] = sort(bdnb1(:,1)*num+bdnb1(:,2));
bdnb2 = bdnb1(ord,:);
ff   = bondnb(fixs,:);
ff1  = floor((ff-1)/nx)+mod(ff-1,nx)*ny+1;
ord1 = zeros(size(ord));
ii0  = find(ismember(bdnb2,ff1,'rows'));
i0 = 1;
for i=1:length(ii0)
    ord1(i0:ii0(i)-i) = i0+i-1:ii0(i)-1;
    i0 = ii0(i)-i+1;
end
ord1(i0:length(ord)-length(ii0))=i0+length(ii0):length(ord);
ord1(length(ord)-length(ii0)+1:length(ord))=ii0;
dirc  = '../Allostery/largerSystems_nov16_beta100/';
nx = 20; ny = 20;
xname  = sprintf('%d',nx);
yname  = sprintf('%d',ny);
num    = ny*nx;
nsp    = round(num*5/2);
nname  = sprintf('%d',nsp);
for r = 1:rp
    rname = sprintf('%d',r);
    dirct = [dirc,'Nx',xname,'_Ny',yname,'_Nsp',nname,'_MCS50000_Tw50000_',rname,'/'];
    filename = [dirct,'configs_Nx',xname,'_Nx',yname,'_Nsp',nname,'_d-1.0e+00_d1-1.0e+00_delta-0.20_kw-1.0e-04_GSL_beta100.00000_MCS50000_tw50000_0.dat'];
    data = dlmread(filename);
    config = zeros(size(data,1),nnb);
    for i=1:size(data,1)
        config(i,ord(ord1(data(i,:)+1)))=1;
    end
    filename = [dirc,cname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_1_',rname,dtype]; 
    dlmwrite(filename,config);
end

%% read the configurations
%
dirc  = '../Allostery/May23/';
nx = 12; ny = 12; rp = 20;
xname  = sprintf('%d',nx);
yname  = sprintf('%d',ny);
num    = ny*nx;
dim    = 2;
nsp    = round(num*5/2);
ew     = zeros(num,1);  % energy stored in weak springs averaged on sites
es     = zeros(num,1);  % energy stored in strong springs averaged on sites
%ewk    = zeros(nwk,1);  % energy stored in weak springs
%estr   = zeros(3*num-2*nx,1);  % energy stored in strong springs
piso   = zeros(num,1);  % probability of the sites in the isostatic network

sigm   = zeros(3*num-2*nx,1);  % probability of spring 
psig   = zeros(num,1);         % prob of spring average for site
nz     = zeros(num,1);  
mp     = zeros(num,1);
pm     = zeros(dim*num,1);
%id = 1:4;
%idd = 1:8;
%nid = 5:num;
nname  = sprintf('%.3f',nsp/(3*num-2*nx)); %
tnam   = sprintf('%.4f',0.4);
nr     = 0;
tic;
for r = 1:rp
    rname = sprintf('%03d',r);
    filename = [dirc,cname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_1_',rname,dtype]; %
    config = dlmread(filename);
    cl    = 0; %floor(size(config,1)/2);
    for tt = 1:size(config,1)  %cl
        stat = config(cl+tt,:);
        sigm = sigm+stat';
        %{
        Qmat   = zeros(dim*num);
        for i=idd
            Qmat(i,i) = 1;
        end
        Mmat = kweak*MmatW+Smatrix'*diag(stat)*Smatrix;
        Qmat(:,nid) = -Mmat(:,nid);
        fext = Mmat*pert;
        disp = pinv(Qmat)*fext; %linsolve(Qmat,fext); % Qmat\fext; %
        displ = zeros(dim*num,1);
        xid = mod(nid,2)==1;
        yid = mod(nid,2)==0;
        dx = mean(disp(nid(xid)));
        dy = mean(disp(nid(yid)));
        %[V,~] = eig(Mmat);
        %
        displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
        %if nsp>2*xl*yl %subtract the first phonon
        %y0 = mean(posy);
        %nm = norm(posy-y0);
        %dx1 = (posy-y0)*displ(1:2:end)/nm;
        %displ(1:2:end) = displ(1:2:end)-dx1*(posy'-y0)/nm;
        %end
        %dispx = displ(1:2:end);
        %dispy = displ(2:2:end);
        mp = mp+(displ(1:2:end).^2+displ(2:2:end).^2)/sum(displ.^2);
        pm = pm+displ/norm(displ);
        %}
        bdn =bondnb(stat>0,:);
        nz   = nz+histc(reshape(bdn,1,[]),1:num)';
        %sites = findIsostatic(it,num-nx+1:num,bdn,[posx;posy]');
        %force = diag(stat)*Smatrix*(pert+disp);
        %fweak = SMatW*(pert+disp);
        %estr  = estr+force.^2/2;
        %ewk   = ewk+fweak.^2/2;
        %piso(sites) = piso(sites)+1;
        %}
        nr    = nr+1;
    end
end
toc;
%estr = estr/nr;
%ewk  = ewk/nr;
piso = piso/nr;
sigm = sigm/nr;
nz   = nz/nr;
pm   = pm/nr;
mp   = mp/nr;
dsig = sigm.*log(sigm/mean(sigm))+(1-sigm).*log((1-sigm)/(1-mean(sigm)));
dsig(sigm==0) = log(1/(1-mean(sigm)));
dsig(sigm==1) = log(1/mean(sigm));
for i=1:num
    id = sum(ismember(bondnb,i),2)>0;
%    es(i) = sum(estr(id))/2;
    psig(i) = mean(dsig(id));
%    id = sum(ismember(bondwk,i),2)>0;
%    ew(i) = sum(ewk(id))/2;
end
%}