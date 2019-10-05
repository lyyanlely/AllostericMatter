%% read the numerical data
dirc  = '../Allostery/May4/';
ename = 'Costs';
cname = 'Config';
dtype = '.dat';
dd    = 1;
dt    = 1;
kweak = 0.0001;
dname = sprintf('%.2f',dd);
tname = sprintf('%.2f',dt);
kname = sprintf('%.4f',kweak);

rp    = 1;
npara = 9;

beta = 1./(0.1:0.2/(npara-1):0.3);

nx  = 12; %
ny  = 12; %[4,6,8,10,12,14,16]; %[3,4,5,6,7,8,9,10,11,12];

lx  = length(ny);

nsp = zeros(lx,1);
nmax= zeros(lx,1);
pos = cell(lx,1);
posb= cell(lx,1);
mz  = cell(lx,1);
vz  = cell(lx,1);

xname = sprintf('%d',nx);
for j = 1:lx
    xl = nx;
    yl = ny(j);
    num    = yl*xl;
    yname  = sprintf('%d',yl);
    nsp(j) = round(xl*yl*4.25/2); %(2*nx*ny(j)-4):(3*nx*ny(j)-2*nx-1);
    nmax(j)= 3*yl*xl-2*xl;
    nname  = sprintf('%d',nsp(j));
    mz{j}  = zeros(npara,num);
    vz{j}  = zeros(npara,num);
    %% construct the embedded network
    bondnb = zeros(3*num,2);
    nnb    = 0;
    posx = zeros(1,num);
    posy = zeros(1,num);
    posbx = zeros(1,3*num);
    posby = zeros(1,3*num);
    for n=1:num
        yy = ceil(n/xl);
        xx = n - (yy-1)*xl;
        posx(n) = mod(xx-1+4,xl)+mod(yy-1,2)/2;
        posy(n) = sqrt(3)/2*yy - sqrt(3)/4;
        %
        if mod(yy,2)
            if mod(xx+mod(floor(yy/2),2),2)
                posx(n) = posx(n);
                posy(n) = posy(n);
            else
                posx(n) = posx(n) + eps*sqrt(3)/2;
                posy(n) = posy(n) + eps/2;
            end
        else
            if mod(xx+mod(yy/2-1,2),2)
                posx(n) = posx(n);
                posy(n) = posy(n) - eps;
            else
                posx(n) = posx(n) - eps*sqrt(3)/2;
                posy(n) = posy(n) + eps/2;
            end
        end
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+1-xl*(xx==xl);  % periodic in horizontal direction
        posbx(nnb) = posx(n)+1/2;
        posby(nnb) = posy(n);
            
        if yy<yl  % the top line is not periodic
            %    if nx<xl
            nnb = nnb+1;
            bondnb(nnb,1) = n;
            bondnb(nnb,2) = n+xl+(1-xl*(xx==xl))*mod(yy-1,2); %-num*(ny==linnum)
            posbx(nnb) = posx(n)+1/4;
            posby(nnb) = posy(n)+sqrt(3)/4;
            %    end
            %    if nx>1
            nnb = nnb+1;
            bondnb(nnb,1) = n;
            bondnb(nnb,2) = n+xl-(1-xl*(xx==1))*mod(yy,2);  %-num*(ny==linnum)
            posbx(nnb) = posx(n)-1/4;
            posby(nnb) = posy(n)+sqrt(3)/4;
            %    end
        end
        
    end
    bondnb = bondnb(1:nnb,:);
    pos{j} = [posx;posy];
    posb{j} = [posbx(1:nnb);posby(1:nnb)];
    %% read different temperature
    for k = 1:npara
        pname  = sprintf('%d',k);
        mztemp = zeros(1,num);
        vztemp = zeros(1,num);
        nr     = 0;
        for r = 1:rp
            rname = sprintf('%03d',r);
            filename = [dirc,cname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',kname,'_',pname,'_',rname,dtype];
            if exist(filename,'file')
                s = dir(filename);
                if s.bytes>1
                    config = dlmread(filename);
                    if min(size(config))<nmax(j)
                        config = reshape(config,[],nmax(j));
                    end
                    cl    = size(config);
                    cl    = cl(1);
                    coord = zeros(cl/2,num);
                    for c = 1:cl/2
                        state = config(c+cl/2,:);
                        [cc,~]=hist(reshape(bondnb(state>0,:),1,[]),1:num);
                        coord(c,:) = cc;
                        if c==1
                            nr       = nr+1;
                        end
                    end
                    mztemp = mztemp+mean(coord);
                    vztemp = vztemp+var(coord);
                end
            end
        end
        mz{j}(k,:) = mztemp/nr;
        vz{j}(k,:) = vztemp/nr;
    end
end