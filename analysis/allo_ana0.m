%% read the numerical data
dirc  = '../Allostery/May4/';
ename = 'Costs';
zname = 'MeanZ';
pname = 'PartRatio';
cname = 'Config';
plnam = 'PrLocal';
zlnam = 'Zlocal';
bfnam = 'Bfactor';
sfnam = 'Sfactor';
dtype = '.dat';
dd    = 1;
dt    = 1;
kweak = 0.0001;
dname = sprintf('%.2f',dd);
tname = sprintf('%.2f',dt);
kname = sprintf('%.4f',kweak);

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
beta = [1e12,1./(0.04:0.015:.16)];

nx  = [4,6,8,10,12];
ny  = [4,6,8,10,12]; %[3,4,5,6,7,8,9,10,11,12];

lx  = length(nx);
npara = length(beta);

num  = zeros(lx,1);
nsp  = cell(lx,1);
cost = cell(lx,1);
hcap = cell(lx,1);
cdz  = cell(lx,1);
part = cell(lx,1);
zloc = cell(lx,1);
ploc = cell(lx,1);
bfac = cell(lx,1);
sfac = cell(lx,1);

for j = 1:lx
    xname  = sprintf('%d',nx(j));
    yname  = sprintf('%d',ny(j));
    num(j) = ny(j)*nx(j);
    nsp{j} = round(num(j)*3/2); %(2*nx*ny(j)-4):(3*nx*ny(j)-2*nx-1);
    ln  = length(nsp{j});
    cost{j} = zeros(ln,npara);
    cdz{j}  = zeros(ln,npara);
    part{j} = zeros(ln,npara);
    zloc{j} = zeros(npara,num(j));
    ploc{j} = zeros(npara,num(j));
    %bfac{j} = zeros(npara,num(j));
    %sfac{j} = zeros(npara,3*num(j)-2*nx(j));
    for k = 1:ln
        nname  = sprintf('%d',nsp{j}(k));
        ctemp  = zeros(1,npara);
        htemp  = zeros(1,npara);
        ztemp  = zeros(1,npara);
        ptemp  = zeros(1,npara);
        zltmp  = zeros(npara,num(j));
        pltmp  = zeros(npara,num(j));
        %bftmp  = zeros(npara,num(j));
        %sftmp  = zeros(npara,3*num(j)-2*nx(j));
        %for n = 1:npara
        %    tnam = sprintf('%.4f',beta(n));
        nr     = 0;
        for r = 1:rp
            rname = sprintf('%03d',r);
            filename = [dirc,ename,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',rname,dtype];  %'_',tnam,
            cordname = [dirc,zname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',rname,dtype];
            partname = [dirc,pname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',rname,dtype];
            zlocname = [dirc,zlnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',rname,dtype];
            prlcname = [dirc,plnam,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',rname,dtype];
            if exist(filename,'file')
                s = dir(filename);
                if s.bytes>1
                    cdata = dlmread(filename);
                    if min(size(cdata))<npara
                        cdata = reshape(cdata,[],npara);
                    else
                        cdata = cdata';
                    end
                    zdata = dlmread(cordname);
                    if min(size(zdata))<npara
                        zdata = reshape(zdata,[],npara);
                    else
                        zdata = zdata';
                    end
                    pdata = dlmread(partname);
                    if min(size(pdata))<npara
                        pdata = reshape(pdata,[],npara);
                    else
                        pdata = pdata';
                    end
                    zldat = dlmread(zlocname);
                    if min(size(zldat))<npara
                        zldat = reshape(zldat,npara,[]);
                    end
                    pldat = dlmread(prlcname);
                    if min(size(pldat))<npara
                        pldat = reshape(pldat,npara,[]);
                    end
                    cl    = size(cdata);
                    cl    = cl(1);
                    for n = 1:npara
                        ctemp(n) = ctemp(n)+mean(cdata(round(cl/2)+1:cl,n));
                        htemp(n) = htemp(n)+var(cdata(round(cl/2)+1:cl,n));
                        ztemp(n) = ztemp(n)+mean(zdata(round(cl/2)+1:cl,n));
                        ptemp(n) = ptemp(n)+mean(pdata(round(cl/2)+1:cl,n));
                        zltmp(n,:) = zltmp(n,:)+zldat(n,:);
                        pltmp(n,:) = pltmp(n,:)+pldat(n,:);
                        %bftmp(n,:) = bftmp(n,:)+bfdat(1,:);
                        %sftmp(n,:) = sftmp(n,:)+sfdat(1,:);
                        if n==1
                            nr       = nr+1;
                        end
                    end
                end
            end
        end
        cost{j}(k,:) = ctemp/nr;
        hcap{j}(k,:) = htemp/nr;
        cdz{j}(k,:)  = ztemp/nr;
        part{j}(k,:) = ptemp/nr;
        zloc{j} = zltmp/nr;
        ploc{j} = pltmp/nr;
        %bfac{j} = bftmp/nr;
        %sfac{j} = sftmp/nr;
    end
    
    %% compute the eigenmodes
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
