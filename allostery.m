function allostery(replica, xl, yl, nsp, kweak)
%replica = 1;
%xl    = 12;
%yl    = 12;
%nsp   = 360;%
dd    = 1;
dt    = 1.;
%kweak = 0.0001;
beta_ini = 0.0001; %[0.01,0.03,0.05];
np    = 4;  % number of sites to perturb
nt    = 4;  % number of sites in target
temprange = length(beta_ini);
iteration = 100000;
paraiter  = 100000;
interv = 1000;
nrec   = round(iteration/interv);
prec   = round(paraiter/interv);
npara  = round(iteration/paraiter);
eps   = 0.2;   % distort size
%confint = 100;
%slow = 10;
warning('off','all');
%correlatetau = 1/exp(1);
%svol  = 1;
%% size of system and size dependent parameters
num   = xl*yl;
dim   = 2;
xbound = xl;
ybound = yl*sqrt(3)/2;
nmax  = 3*num-2*xl;
%pp    = nsp/nmax;   %fraction of occupied springs

%tol   = 1/65536/num;
%maxit = 16*num;

%beta  = zeros(1,temprange);
%temprange = length(beta);
nflip = ones(1,temprange);
nf    = ones(1,temprange);
%{
if temprange>3
    b0   = beta_ini^(1/(temprange-2));
    if b0>1
        beta(2:temprange-1) = b0.^(0:temprange-3);
    else
        beta(2:temprange-1) = b0.^(3-temprange:0);
    end
    beta(temprange) = 1e12;
elseif temprange==3
    beta = [1,beta_ini,1e12];
elseif temprange==2
    beta = [beta_ini,1e12];
else
    beta = 1e12;
end
%}
beta  = 1./beta_ini;
for tt = 1:temprange
    nflip(tt) = temprange+1-tt;
end

%% define perturb and target displacements
id  = (1:np);    %floor((xl-np)/2)+ indices of the displaced particles
idd = dim*id(1)-1:dim*id(end);  %
it  = (yl-1)*xl+(1:nt);  %floor((xl-nt)/2)+ indices of the target particles
itt = dim*it(1)-1:dim*it(end);
nid = setxor(1:dim*num,idd);
nit = setxor(1:dim*num,itt);
nidt= setxor(1:dim*num,[idd,itt]);
pert = zeros(dim*num,1);  % perturb displacement
targ = zeros(dim*num,1);  % target displacement
fixs = zeros(1,(np+nt)*2);
cf   = 0;
%% construct the embedded network
bondnb = zeros(3*num,2);
nnb    = 0;
bondwk = zeros(6*num,2);  % to all next neareat neighbors
nwk    = 0;
posx = zeros(1,num);
posy = zeros(1,num);
for n=1:num
    ny = ceil(n/xl);
    nx = n - (ny-1)*xl;
    posx(n) = nx-1+mod(ny-1,2)/2;
    posy(n) = sqrt(3)/2*ny - sqrt(3)/4;
    %
    if mod(ny,2)
        if mod(nx+mod(floor(ny/2),2),2)
            posx(n) = posx(n);
            posy(n) = posy(n);
        else
            posx(n) = posx(n) + eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    else
        if mod(nx+mod(ny/2-1,2),2)
            posx(n) = posx(n);
            posy(n) = posy(n) - eps;
        else
            posx(n) = posx(n) - eps*sqrt(3)/2;
            posy(n) = posy(n) + eps/2;
        end
    end
    %}
    %theta = rand(1);
    %posx(n) = posx(n) + eps*cos(theta);
    %posy(n) = posy(n) + eps*sin(theta);
    if nx<xl
    nnb = nnb+1;
    bondnb(nnb,1) = n;
    bondnb(nnb,2) = n+1-xl*(nx==xl);  % periodic in horizontal direction
    end
    % fix the connections between the perturbing nodes and target nodes
    if ismember(n,id)&&ismember(n+1-xl*(nx==xl),id)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ismember(n,it)&&ismember(n+1-xl*(nx==xl),it)
        cf = cf+1;
        fixs(cf) = nnb;
    end
    if ny<yl  % the top line is not periodic
        if nx<xl||mod(ny,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl+(1-xl*(nx==xl))*mod(ny-1,2); %-num*(ny==linnum)
        end
        if nx>1||mod(ny+1,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl-(1-xl*(nx==1))*mod(ny,2);  %-num*(ny==linnum)
        end
    end
    % weak bonds connect to the next nearest and next next nearest
    % next nearest neighbors
    if ny<yl
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl+1-xl*(nx==xl)+(1-xl*(nx==xl-1))*mod(ny-1,2);  %-num*(ny==linnum)
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+xl-1+xl*(nx==1)-(1-xl*(nx==2))*mod(ny,2);  %-num*(ny==linnum)
    end
    if ny<yl-1
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl;
        % next-next nearest neighbors
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl+1-xl*(nx==xl);
        nwk = nwk+1;
        bondwk(nwk,1) = n;
        bondwk(nwk,2) = n+2*xl-1+xl*(nx==1);
    end
    nwk = nwk+1;
    bondwk(nwk,1) = n;
    bondwk(nwk,2) = n+2-xl*(nx>=xl-1);
end
bondnb = bondnb(1:nnb,:);
bondwk = bondwk(1:nwk,:);
fixs   = fixs(1:cf);
pert(2*id(1)) = -1*dd;
targ(2*it(1)) = 1*dt;
for i = 2:np
    dx  = posx(id(i))-posx(id(i-1));
    dy  = posy(id(i))-posy(id(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    pert(2*id(i)-1) = pert(2*id(i-1)-1)*(dx^2-dy^2)+2*pert(2*id(i-1))*dx*dy;
    pert(2*id(i))   = pert(2*id(i-1))*(dy^2-dx^2)+2*pert(2*id(i-1)-1)*dx*dy;
end
for i = 2:nt
    dx  = posx(it(i))-posx(it(i-1));
    dy  = posy(it(i))-posy(it(i-1));
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    targ(2*it(i)-1) = targ(2*it(i-1)-1)*(dx^2-dy^2)+2*targ(2*it(i-1))*dx*dy;
    targ(2*it(i))   = targ(2*it(i-1))*(dy^2-dx^2)+2*targ(2*it(i-1)-1)*dx*dy;
end

% perturb of one bond on boundary
ptmp = zeros(xl,dim*num);
id0  = zeros(xl,2);
idd0 = zeros(xl,dim*2);
for i=1:xl
    i0 = i;
    id0(i,:) = [i0,mod(i0,xl)+1]; %
    idd0(i,:)= [2*id0(i,:)-1,2*id0(i,:)];
    ptmp(i,2*i0-1) = 0;
    ptmp(i,2*i0)   = 1;
    dx  = posx(id0(i,2))-posx(id0(i,1));
    dy  = posy(id0(i,2))-posy(id0(i,1));
    if dx>xl/2
        dx = dx-xl;
    elseif dx<-xl/2
        dx = dx+xl;
    end
    dr  = sqrt(dx^2+dy^2);
    dx  = dx/dr;
    dy  = dy/dr;
    ptmp(i,idd0(i,2)) = ptmp(i,idd0(i,1))*(dx^2-dy^2)+2*ptmp(i,idd0(i,3))*dx*dy;
    ptmp(i,idd0(i,4)) = ptmp(i,idd0(i,3))*(dy^2-dx^2)+2*ptmp(i,idd0(i,1))*dx*dy;
end
%% Compute Mmatrix
bpx  = zeros(1,nnb);
bpy  = zeros(1,nnb);
delx = zeros(1,nnb);
dely = zeros(1,nnb);
delr = zeros(1,nnb);
idx = zeros(1,4*nnb);
jdx = zeros(1,4*nnb);
val = zeros(1,4*nnb);
for n = 1:nnb
    ii = bondnb(n,1);
    jj = bondnb(n,2);
    delx(n) = posx(ii) - posx(jj);
    dely(n) = posy(ii) - posy(jj);
    if delx(n)>xbound/2   % periodic
        delx(n) = delx(n) - xbound;
    elseif delx(n)<-xbound/2
        delx(n) = delx(n) + xbound;
    end
    bpx(n) = posx(jj)+delx(n)/2;
    bpy(n) = posy(jj)+dely(n)/2;
    %if dely(n)>ybound/2  % not periodic
    %    dely(n) = dely(n) - ybound;
    %elseif dely(n)<-ybound/2
    %    dely(n) = dely(n) + ybound;
    %end
    delr(n) = sqrt(delx(n)^2+dely(n)^2);
    idx(4*(n-1)+1) = n;
    jdx(4*(n-1)+1) = 2*(ii-1)+1;
    val(4*(n-1)+1) = delx(n)/delr(n);
    idx(4*(n-1)+2) = n;
    jdx(4*(n-1)+2) = 2*(ii-1)+2;
    val(4*(n-1)+2) = dely(n)/delr(n);
    idx(4*(n-1)+3) = n;
    jdx(4*(n-1)+3) = 2*(jj-1)+1;
    val(4*(n-1)+3) = -delx(n)/delr(n);
    idx(4*(n-1)+4) = n;
    jdx(4*(n-1)+4) = 2*(jj-1)+2;
    val(4*(n-1)+4) = -dely(n)/delr(n);
end
Smatrix = sparse(idx,jdx,val,nnb,dim*num);
%Tmatrix = Smatrix';
% of the weak network
bpxw = zeros(1,nnb);
bpyw = zeros(1,nnb);
delxw = zeros(1,nwk);
delyw = zeros(1,nwk);
delrw = zeros(1,nwk);
idxw = zeros(1,4*nwk);
jdxw = zeros(1,4*nwk);
valw = zeros(1,4*nwk);
for n = 1:nwk
    ii = bondwk(n,1);
    jj = bondwk(n,2);
    delxw(n) = posx(ii) - posx(jj);
    delyw(n) = posy(ii) - posy(jj);
    if delxw(n)>xbound/2
        delxw(n) = delxw(n) - xbound;
    elseif delxw(n)<-xbound/2
        delxw(n) = delxw(n) + xbound;
    end
    %if delyw(n)>ybound/2
    %    delyw(n) = delyw(n) - ybound;
    %elseif delyw(n)<-ybound/2
    %    delyw(n) = delyw(n) + ybound;
    %end
    bpxw(n) = posx(ii)+delxw(n)/2;
    bpyw(n) = posy(ii)+delyw(n)/2;
    delrw(n) = sqrt(delxw(n)^2+delyw(n)^2);
    idxw(4*(n-1)+1) = n;
    jdxw(4*(n-1)+1) = 2*(ii-1)+1;
    valw(4*(n-1)+1) = delxw(n)/delrw(n);
    idxw(4*(n-1)+2) = n;
    jdxw(4*(n-1)+2) = 2*(ii-1)+2;
    valw(4*(n-1)+2) = delyw(n)/delrw(n);
    idxw(4*(n-1)+3) = n;
    jdxw(4*(n-1)+3) = 2*(jj-1)+1;
    valw(4*(n-1)+3) = -delxw(n)/delrw(n);
    idxw(4*(n-1)+4) = n;
    jdxw(4*(n-1)+4) = 2*(jj-1)+2;
    valw(4*(n-1)+4) = -delyw(n)/delrw(n);
end
SMatW = sparse(idxw,jdxw,valw,nwk,dim*num);
MmatW = SMatW'*SMatW;


%{
vecShr = zeros(1, 3*num);
for i=1:3*num
    vecShr(i) = (Smatrix(i,2*bondnb(i)-1)*Smatrix(i,2*bondnb(i)));
end

hlayer = zeros(linnum,2*linnum);
vlayer = zeros(linnum,2*linnum);
for i=1:linnum
    for j=1:linnum
        hlayer(i, 2*(j-1)+1) = 3*linnum*(i-1)+3*(j-1)+2;
        hlayer(i, 2*j)       = 3*linnum*(i-1)+3*j;
        vlayer(i, 2*(j-1)+1) = 3*linnum*(j-1)+3*(i-1)+1;
        vlayer(i, 2*j)       = 3*linnum*(j-1)+3*(i-1)+2+4*mod(j,2)-mod(j,2)*floor(i/linnum)*3*linnum;
    end
end
%}
%% data structure
s = rng('shuffle');
i0    = 2;
nz     = zeros(temprange,dim*num);
ploc   = zeros(temprange,num);    % local participation ratio
pplc   = zeros(temprange,num*dim);    % local participation ratio
zloc   = zeros(temprange,num);    % local coordination number
estr   = zeros(temprange,nnb);    % energy in strong springs
ewk    = zeros(temprange,nwk);    % energy in weak springs
bfact  = zeros(temprange,num);    % b-factor
sfact  = zeros(temprange,nnb);    % s-factor
cloc   = zeros(temprange,nnb);    % cost from changing the state
nc     = zeros(temprange,1);      % count the number of contributions
cdata  = zeros(temprange,nrec);   % cost
pdata  = zeros(temprange,nrec);   % partition ratio
zdata  = zeros(temprange,nrec);   % coordination of the most participated 
estim  = zeros(temprange,nrec);   % energy cost due to stimulus
etarg  = zeros(temprange,nrec);   % energy cost due to target
esttg  = zeros(temprange,nrec);   % energy cost due to cooperativity of stimulus and target
config = zeros(nrec,nnb,temprange);
dd     = zeros(20,nrec);
qq     = zeros(20,nrec);
dd1    = zeros(20,nrec);
qq1    = zeros(20,nrec);
resps  = zeros(xl,num);
respv  = zeros(xl,num*dim);
for tt = 1:temprange
    nf(tt) = (min([nflip(tt),nsp,nnb-nsp]));
end

%% Initialize the configurations
cost   = zeros(temprange,1);
prat   = zeros(temprange,1);
coord  = zeros(temprange,1);
occupy = zeros(temprange,nsp-cf);
vacant = zeros(temprange,nnb-nsp);
state  = zeros(temprange,nnb);
%Mmat   = zeros(dim*num);
Qmat   = zeros(dim*num);
for tt = 1:temprange
    occupy(tt,:) = randsample(setdiff(1:nnb,fixs), nsp-cf)'; %
    vacant(tt,:) = setxor(setdiff(1:nnb,fixs), occupy(tt,:)); %
    state(tt,union(occupy(tt,:),fixs)) = 1; %
    Mmat = kweak*MmatW+Smatrix'*diag(state(tt,:))*Smatrix;
    for i=idd
        Qmat(i,i) = 1;
    end
    Qmat(:,nid) = -Mmat(:,nid);
    fext = Mmat*pert;
    disp = linsolve(Qmat,fext);  %Qmat\fext; % 
    if any(isnan(disp)|abs(disp)>1e5)
        disp = pinv(Qmat)*fext;
    end
    diff = disp(itt)-targ(itt);
    cost(tt) = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
    %compute the coordination for each node
    [cc,~]=hist(reshape(bondnb(union(fixs,occupy(tt,:)),:),1,[]),1:num);
    for d = 1:dim
        nz(tt,d:dim:end) = cc;
    end
    %
    displ = zeros(dim*num,1);
    xid = mod(nid,2)==1;
    yid = mod(nid,2)==0;
    dx = mean(disp(nid(xid)));
    dy = mean(disp(nid(yid)));
    displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
    if nsp>2*xl*yl %subtract the first phonon
        y0 = mean(posy);
        nm = norm(posy-y0);
        dx1 = (posy-y0)*displ(1:2:end)/nm;
        displ(1:2:end) = displ(1:2:end)-dx1*(posy'-y0)/nm;
    end
    prat(tt) = sum(displ(nz(tt,:)>2).^2)^2/num/dim/sum(displ(nz(tt,:)>2).^4);
    coord(tt)= nz(tt,nz(tt,:)>2)*(displ(nz(tt,:)>2).^2)/sum(displ(nz(tt,:)>2).^2);
    %}
end
%{
pebble = zeros(2,num);
block  = -2*ones(1,num);
btree  = zeros(1,num);
shell  = zeros(1,num);
ndep   = 0;
tag    = 0;
coverindex = zeros(1,numcon);  %a number labels a pebble from which site covers, 0 no covering pebble overconstraint bond
for n=1:numcon
    newbd = bondnb(occupy(n),:);
    [ndep_temp,tag,pebble,shell,block,btree] = check(newbd(1),newbd(2),ndep,tag,pebble,shell,block,btree);
    if ndep_temp>ndep
        coverindex(n) = 1;
        ndep = ndep_temp;
    end
end
Nr_temp = ndep;
%}
%% iteration of Monte Carlo
tic;
for t = 1:npara
    for tt = 1:temprange
        %ptarg = targ;
        pbeta = beta(tt);
        pnf   = randi(nf(tt));
        pocc  = occupy(tt,:);
        pvac  = vacant(tt,:);
        pstat = state(tt,:);
        pcost = cost(tt);
        pr    = prat(tt);
        pz    = coord(tt);
        pcdat = zeros(1,prec);
        partr = zeros(1,prec);
        pcoord= zeros(1,prec);
        es    = zeros(1,prec);
        et    = zeros(1,prec);
        est   = zeros(1,prec);
        %partr(1) = prat(tt);
        pconf = zeros(prec,nnb);
        pQmat = zeros(dim*num);
        for i=idd
            pQmat(i,i) = 1;
        end
        for ii = 1:paraiter
            ctemp = pcost;
            statetemp = pstat;
            % new state
            old = randsample(pocc,pnf);
            if length(pvac)>1
                new = randsample(pvac,pnf);
            else
                new = pvac;
            end
            pstat(old) = 0;
            pstat(new) = 1;
            pMmat = kweak*MmatW+Smatrix'*diag(pstat)*Smatrix;
            pQmat(:,nid) = -pMmat(:,nid);
            fext = pMmat*pert;
            disp = linsolve(pQmat,fext);  %pQmat\fext; %
            if any(isnan(disp)|abs(disp)>1e5)
                disp = pinv(pQmat)*fext; %
            end
            diff = disp(itt)-targ(itt);
            pcost = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
            % metropolis criterion
            if rand(1)<exp(-pbeta*(pcost-ctemp))
                pocc(ismember(pocc,old)) = new;
                pvac(ismember(pvac,new)) = old;
                displ = zeros(dim*num,1);
                xid = mod(nid,2)==1;
                yid = mod(nid,2)==0;
                dx = mean(disp(nid(xid)));
                dy = mean(disp(nid(yid)));
                displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
                if nsp>2*xl*yl %subtract the first phonon
                    y0 = mean(posy);
                    nm = norm(posy-y0);
                    dx1 = (posy-y0)*displ(1:2:end)/nm;
                    displ(1:2:end) = displ(1:2:end)-dx1*(posy'-y0)/nm;
                end
                pr  = sum(displ(nz(tt,:)>2).^2)^2/num/dim/sum(displ(nz(tt,:)>2).^4);
                remv = reshape(bondnb(old,:),1,[]);
                badd = reshape(bondnb(new,:),1,[]);
                for d=1:dim
                    for i=remv
                        nz(tt,dim*(i-1)+d) = nz(tt,dim*(i-1)+d)-1;
                    end
                    for i=badd
                        nz(tt,dim*(i-1)+d) = nz(tt,dim*(i-1)+d)+1;
                    end
                end
                pz  = nz(tt,nz(tt,:)>2)*(displ(nz(tt,:)>2).^2)/sum(displ(nz(tt,:)>2).^2);
            else
                pcost = ctemp;
                pstat = statetemp;
            end
            if mod(ii,10000)==0
                1;
            end
            % record
            if mod(ii,interv)==0
                pcdat(ii/interv)   = pcost;
                pconf(ii/interv,:) = pstat;
                partr(ii/interv)   = pr;
                pcoord(ii/interv)  = pz;
                %
                disp(idd) = pert(idd);
                force = diag(pstat)*Smatrix*(pert+disp);
                fweak = SMatW*(pert+disp);
                es(ii/interv)      = 0.5*disp'*pMmat*disp;
                targp  = zeros(size(targ));
                targp(itt(1:2:end)) = targ(itt(1:2:end))+mean(disp(itt(1:2:end))-targ(itt(1:2:end)));
                targp(itt(2:2:end)) = targ(itt(2:2:end))+mean(disp(itt(2:2:end))-targ(itt(2:2:end)));
                Qmat   = zeros(dim*num);
                for i=itt
                    Qmat(i,i) = 1;
                end
                Qmat(:,nit) = -pMmat(:,nit);
                fext = pMmat*targp;
                disp = linsolve(Qmat,fext);  %pQmat\fext; %
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext; %
                end
                disp(itt) = targp(itt);
                et(ii/interv) = 0.5*disp'*pMmat*disp;
                Qmat   = zeros(dim*num);
                for i=itt
                    Qmat(i,i) = 1;
                end
                for i=idd
                    Qmat(i,i) = 1;
                end
                Qmat(:,nidt) = -pMmat(:,nidt);
                fext = pMmat*(pert+targp);
                disp = linsolve(Qmat,fext);  %pQmat\fext; %
                if any(isnan(disp)|abs(disp)>1e5)
                    disp = pinv(Qmat)*fext; %
                end
                disp(idd) = pert(idd);
                disp(itt) = targp(itt);
                est(ii/interv) = 0.5*disp'*pMmat*disp;
                %}
                % compute the eigenvalues
                pMmat = kweak*MmatW+Smatrix'*diag(pstat)*Smatrix;
                Mmat  = Smatrix'*diag(pstat)*Smatrix;
                [V,D] = eig(full(pMmat));
                [Vr,Dr] = eig(full(Mmat));
                d = diag(D);
                d1= diag(Dr);
                dd(:,ii/interv) = sqrt(d(i0+(1:20)));
                qq(:,ii/interv) = abs(V(:,i0+(1:20))'*displ)/norm(displ);
                dd1(:,ii/interv) = sqrt(d1(i0+(1:20)));
                qq1(:,ii/interv) = abs(Vr(:,i0+(1:20))'*displ)/norm(displ);
                % compute the perturbation
                %{
                pnbc = zeros(1,nnb);
                for i=1:nnb
                    ppstat = pstat;
                    ppstat(i) = mod(pstat(i)+1,2);
                    ppMmat = kweak*MmatW+Smatrix'*diag(ppstat)*Smatrix;
                    pQmat(:,nid) = -ppMmat(:,nid);
                    fext = ppMmat*pert;
                    disp = linsolve(pQmat,fext);  %pQmat\fext; %
                    if any(isnan(disp)|abs(disp)>1e5)
                        disp = pinv(pQmat)*fext; %
                    end
                    diff = disp(itt)-targ(itt);
                    pnbc(i) = sqrt(sum(diff.^2)-sum(diff.*circshift(diff,dim)));
                end
                %}
                if t>npara/2
                    nc(tt) = nc(tt)+1;
                    zloc(tt,:) = zloc(tt,:)+nz(tt,1:2:end);
                    ddp = displ(nz(tt,:)>2)';
                    ploc(tt,nz(tt,1:2:end)>2) = ploc(tt,nz(tt,1:2:end)>2)+(ddp(1:2:end).^2+ddp(2:2:end).^2)/sum(ddp.^2);
                    pplc(tt,nz(tt,:)>2) = pplc(tt,nz(tt,:)>2)+ddp/norm(ddp);
                    %cloc(tt,:) = cloc(tt,:)+pnbc;
                    estr(tt,:) = estr(tt,:)+force'.^2/2;
                    ewk(tt,:)  = ewk(tt,:)+fweak'.^2/2;
                    %{
                    iMmat = pinv(pMmat);
                    bu  = diag(iMmat)';
                    bfact(tt,:) = bfact(tt,:)+bu(1:2:end)+bu(2:2:end);
                    for b = 1:nnb
                        bij = bondnb(b,:);
                        sfact(tt,b) = sfact(tt,b)+bu(2*bij(1)-1)+bu(2*bij(1))+bu(2*bij(2)-1)+bu(2*bij(2))...
                            -iMmat(2*bij(1)-1,2*bij(2)-1)-iMmat(2*bij(2)-1,2*bij(1)-1)-iMmat(2*bij(1),2*bij(2))-iMmat(2*bij(2),2*bij(1));
                    end
                    % compute the different response
                    for i=1:xl
                        Qmat = zeros(dim*num);
                        for j=idd0(i,:)
                            Qmat(j,j) = 1;
                        end
                        Qmat(:,setdiff(1:dim*num,idd0(i,:))) = -pMmat(:,setdiff(1:dim*num,idd0(i,:)));
                        fext = pMmat*ptmp(i,:)';
                        disp = linsolve(Qmat,fext);  %pQmat\fext; %
                        if any(isnan(disp)|abs(disp)>1e5)
                            disp = pinv(Qmat)*fext; %
                        end
                        disp(idd0(i,:)) = ptmp(i,idd0(i,:))';
                        displ = zeros(dim*num,1);
                        dx = mean(disp(1:2:end));
                        dy = mean(disp(2:2:end));
                        displ(1:2:end) = disp(1:2:end)-dx;
                        displ(2:2:end) = disp(2:2:end)-dy;
                        ddp = displ(nz(tt,:)>2)';
                        resps(i,nz(tt,1:2:end)>2) = resps(i,nz(tt,1:2:end)>2)+(ddp(1:2:end).^2+ddp(2:2:end).^2)/sum(ddp.^2);
                        respv(i,nz(tt,:)>2)       = respv(i,nz(tt,:)>2)+ddp/norm(ddp);
                    end
                    %}
                end
            end
        end
        occupy(tt,:) = pocc;
        vacant(tt,:) = pvac;
        state(tt,:)  = pstat;
        cost(tt)     = pcost;
        prat(tt)     = pr;
        coord(tt)    = pz;
        cdata(tt,(t-1)*prec+(1:prec)) = pcdat;
        config((t-1)*prec+(1:prec),:,tt) = pconf;
        pdata(tt,(t-1)*prec+(1:prec)) = partr;
        zdata(tt,(t-1)*prec+(1:prec)) = pcoord;
        estim(tt,(t-1)*prec+(1:prec)) = es;
        etarg(tt,(t-1)*prec+(1:prec)) = et;
        esttg(tt,(t-1)*prec+(1:prec)) = est;
    end
    % switch paratemp data
    if temprange>1
        for tt = 1:temprange-1
            for t1 = 1:temprange-tt
                if rand(1)<exp((beta(t1)-beta(t1+1))*(cost(t1)-cost(t1+1)))
                    occtemp = occupy(t1,:);
                    vactemp = vacant(t1,:);
                    statemp = state(t1,:);
                    costemp = cost(t1);
                    prtemp  = prat(t1);
                    nztemp  = nz(t1,:);
                    ztemp   = coord(t1);
                    occupy(t1,:) = occupy(t1+1,:);
                    vacant(t1,:) = vacant(t1+1,:);
                    state(t1,:)  = state(t1+1,:);
                    cost(t1)     = cost(t1+1);
                    prat(t1)     = prat(t1+1);
                    nz(t1,:)     = nz(t1+1,:);
                    coord(t1)    = coord(t1+1);
                    occupy(t1+1,:) = occtemp;
                    vacant(t1+1,:) = vactemp;
                    state(t1+1,:)  = statemp;
                    cost(t1+1)     = costemp;
                    prat(t1+1)     = prtemp;
                    nz(t1+1,:)     = nztemp;
                    coord(t1+1)    = ztemp;
                end
            end
        end
    end
end
toc;
for tt=1:temprange
    zloc(tt,:) = zloc(tt,:)/nc(tt);
    ploc(tt,:) = ploc(tt,:)/nc(tt);
    pplc(tt,:) = pplc(tt,:)/nc(tt);
    cloc(tt,:) = cloc(tt,:)/nc(tt);
    estr(tt,:) = estr(tt,:)/nc(tt);
    ewk(tt,:)  = ewk(tt,:)/nc(tt);
    bfact(tt,:) = bfact(tt,:)/nc(tt);
    sfact(tt,:) = sfact(tt,:)/nc(tt);
end
%resps = resps/nc;
%respv = respv/nc;

dirc  = './';
cname = 'Costs';
%sname = 'Response';
%svnam = 'ResponseV';
pname = 'PartRatio';
zname = 'MeanZ';
ppnam = 'PrLocal';
dpnam = 'DispMap';
zznam = 'Zlocal';
bfnam = 'Bfactor';
sfnam = 'Sfactor';
ncnam = 'ChangeCost';
esnam = 'Stimulus';
etnam = 'Target';
stnam = 'Cooperativity';
estnm = 'Estrong';
ewknm = 'Eweak';
xname = sprintf('%d', xl);
yname = sprintf('%d', yl);
nsname = sprintf('%d',nsp);
dname = sprintf('%.2f',dd);
tname = sprintf('%.2f',dt);
wfname = sprintf('%.4f', kweak);
bname  = sprintf('%.4f',beta_ini);
exname = sprintf('%03d', replica);
costname = [cname '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(costname, cdata, '\t');
prname = [pname '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(prname, pdata, '\t');
mzname = [zname '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(mzname, zdata, '\t');
esname = [esnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(esname, estim, '\t');
etname = [etnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(etname, etarg, '\t');
stname = [stnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(stname, esttg, '\t');
estnme = [estnm '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(estnme, estr, '\t');
ewknme = [ewknm '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(ewknme, ewk, '\t');
plname = [ppnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(plname, ploc, '\t');
dpname = [dpnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(dpname, pplc, '\t');
zlname = [zznam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(zlname, zloc, '\t');
bfname = [bfnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(bfname, bfact, '\t');
sfname = [sfnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(sfname, sfact, '\t');
sfname = [ncnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
dlmwrite(sfname, cloc, '\t');
%cfname = [sname '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
%dlmwrite(cfname, resps, '\t');
%cfname = [svnam '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' exname '.dat'];
%dlmwrite(cfname, respv, '\t');
%for tt = 1:temprange
%    ttname = sprintf('%d',tt);
%    confname = [dirc,sname '_' xname '_' yname '_' nsname '_' dname '_' tname '_' bname '_' wfname '_' ttname '_' exname '.dat'];
%    dlmwrite(confname, config(:,:,tt));
%end
