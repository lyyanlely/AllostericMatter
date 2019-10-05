%function ReadTriangular(nx,ny,numcon)
nx = 8;
ny = 12;
num= nx*ny;
nmax = 10000;
numcon = 240;

nm = 3*num-2*nx;
nc = nchoosek(nm,numcon);
%if nc>nmax
%    ncomb = nmax;
%    samp  = randsample(1:nc,ncomb);
%else
    ncomb = nc;
    samp  = 1:nc;
%end

pbstr = zeros(1,ncomb);
rperc = zeros(1,ncomb);
pbsp  = zeros(1,ncomb);
perco = zeros(1,ncomb);
r  = 40586;

dirc  = './'; %/Allostery/generateNetworks
nname = 'network';
xname = sprintf('Nx%d',nx);
yname = sprintf('Ny%d',ny);
sname = sprintf('Nsp%d',numcon);
dtype = '.dat';

%for r = 1:ncomb
    replica = samp(r);
rname = sprintf('n%d',replica-1);
filename = [dirc,nname,'_',xname,'_',yname,'_',sname,'_',rname,dtype];
data  = dlmread(filename);

% extract data
pos = data(:,2:3);
nb = zeros(numcon,2);
cc = 0;
for i = 1:num
    for ii = 4:9
        j = data(i,ii)+1;
        if j>i
            cc= cc+1;
            nb(cc,:) = [i,j];
        end
    end
end
if cc~=numcon
    error('wrong!');
end

xlay   = zeros(ny,nx);
%ylay   = zeros(2*nx,ceil(ny/2));
for i=1:ny
    xlay(i,:) = find(pos(:,2)>sqrt(3)/2*(i-1.5)&pos(:,2)<sqrt(3)/2*(i-0.5))';
end
%for i=1:2*nx
%    ylay(i,:) = find(pos(:,1)>0.5*(i-1.5)&pos(:,1)<0.5*(i-0.5))';
%end
%hlayer = zeros(2*nx,2*ny,2);
vlayer = zeros(2*nx,2,ny-1);
for i=1:nx
    for j=1:ny-1
        %hlayer(i, 2*(j-1)+1) = 3*linnum*(i-1)+3*(j-1)+2;
        %hlayer(i, 2*j)       = 3*linnum*(i-1)+3*j;
        i1 = xlay(j,i);
        j1 = find(min(abs([(pos(:,1)-pos(i1,1))';(pos(:,1)-pos(i1,1)+nx)';(pos(:,1)-pos(i1,1)-nx)']+0.5))<0.5&...
            min(abs([(pos(:,2)-pos(i1,2))';(pos(:,2)-pos(i1,2)+ny*sqrt(3)/2)';(pos(:,2)-pos(i1,2)-ny*sqrt(3)/2)']-sqrt(3)/2))<0.5*sqrt(3)/2);
        j2 = find(min(abs([(pos(:,1)-pos(i1,1))';(pos(:,1)-pos(i1,1)+nx)';(pos(:,1)-pos(i1,1)-nx)']-0.5))<0.5&...
            min(abs([(pos(:,2)-pos(i1,2))';(pos(:,2)-pos(i1,2)+ny*sqrt(3)/2)';(pos(:,2)-pos(i1,2)-ny*sqrt(3)/2)']-sqrt(3)/2))<0.5*sqrt(3)/2);
        vlayer(2*(i-1)+1,:,j) = [min(i1,j1),max(i1,j1)];
        vlayer(2*i,:,j)       = [min(i1,j2),max(i1,j2)];
    end
end

% pebble game
% find the stressed cluster
pebble = zeros(2,num);
block  = -2*ones(1,num);
btree  = zeros(1,num);
shell  = zeros(1,num);
stresd = zeros(1,numcon);
ndep   = 0;
tag    = 0;
coverindex = zeros(1,numcon); % if a bond is not covered by pebble, dependent, the coverindex is 1 otherwise is 0.
for n=1:numcon
    newbd = nb(n,:);
    jmax = 0;
    [ndep_temp,tag,pebble,shell,block,btree,jmax] = check(newbd(1),newbd(2),ndep,tag,pebble,shell,block,btree,jmax);
    if ndep_temp>ndep
        pbpath = shell(1:jmax);
        
        ocindx = setdiff(stresd,0);
        if isempty(ocindx)==0
            for tags=ocindx
                mbds  = find(stresd==tags); %marked bonds
                msts  = union(reshape(nb(mbds,:),1,[]),[]); %marked sites
                if length(intersect(msts,pbpath))>2
                    pbpath = union(pbpath,msts);
                end
            end
        end
        bonds  = all(ismember(nb,pbpath),2); %ismember(occupy,find(ismember(ismember(bondnb,pbpath),[1 1],'rows')));
        stresd(bonds) = tag;
        coverindex(n) = 1;
        ndep = ndep_temp;
    end
end
ocindx = unique(setdiff(stresd,0));

%spnx  = zeros(1,numcon);
spny  = zeros(1,numcon);
span  = zeros(1,numcon);
n = 1;
for tags=ocindx
    rclust = nb((stresd==tags),:);
    %if length(rclust)>=nx
    %    spnx(n) = any(ismember(rclust,hlayer(1,:,:),'rows'));
    %    for ii=2:2*nx
    %        spnx(n) = spnx(n)*any(ismember(rclust,hlayer(ii,:,:),'rows'));
    %    end
    %end
    if length(rclust)>=ny-1
        spny(n) = any(ismember(rclust,vlayer(:,:,1),'rows'));
        for ii=2:ny-1
            spny(n) = spny(n)*any(ismember(rclust,vlayer(:,:,ii),'rows'));
        end
    end
    span(n)  = spny(n); %spnx(n)*
    n = n + 1;
end
strdinsp = [];
n = 1;
for tags=ocindx
    if span(n)>0
        strdinsp = union(strdinsp,find(stresd==tags));
    end
    n = n + 1;
end
sbinsp = length(strdinsp);
pbstr(r) = sbinsp/numcon;
rperc(r) = sign(sum(span));

% identify the rigid clusters, mark gives the cluster label of each bond,
% nrigid is the number of bonds in cluster, rigidset records the labels of bonds
mark     = zeros(1,numcon);
nrigid   = zeros(1,numcon);
rigidset = zeros(numcon);
[~,nrigid,rigidset] = idRigid(mark,nrigid,rigidset,numcon,nb,tag,pebble,shell,block,btree);

% spanning: a path appears in every layer horizontally or vertically.
%spnx  = zeros(1,numcon);
spny  = zeros(1,numcon);
span  = zeros(1,numcon);
n = 1;
for n = 1:numcon
    rclust = nb(rigidset(n,1:nrigid(n)),:);
    %if nrigid(n)>=nx
    %    spnx(n) = any(ismember(rclust,hlayer(1,:,:),'rows'));
    %    for ii=2:2*nx
    %        spnx(n) = spnx(n)*any(ismember(rclust,hlayer(ii,:,:),'rows'));
    %    end
    %end
    if nrigid(n)>=ny-1
        spny(n) = any(ismember(rclust,vlayer(:,:,1),'rows'));
        for ii=2:ny-1
            spny(n) = spny(n)*any(ismember(rclust,vlayer(:,:,ii),'rows'));
        end
    end
    span(n)  = spny(n); %spnx(n)*
end
bondinsp = [];
sitenum  = zeros(1,numcon);
%MatStr   = Smatrix'*diag(state_temp)*Smatrix;
for n=1:numcon
    if span(n)==1
        bondinsp = union(bondinsp,rigidset(n,1:nrigid(n)));
    end
    sitenum(n) = length(union(reshape(nb(rigidset(n,1:nrigid(n)),:),1,[]),[]));
end
nbinsp = length(bondinsp);
pbsp(r) = nbinsp/numcon;
perco(r)= sign(sum(span));
%end

dirc  = '../Allostery/';
pname = 'Probability';
filename = [dirc,pname,'_',xname,'_',yname,'_',sname,dtype];
dlmwrite(filename,[mean(pbstr),mean(rperc),mean(pbsp),mean(perco)]);
%
% plot
ccode1 = [210/255 200/255 100/255];
ccode2 = ccode1*0.5;
axis equal off
hold all
for c = 1:cc
    i = nb(c,1);
    j = nb(c,2);
    dx = pos(j,1)-pos(i,1);
    if ismember(c,find(stresd>0))
        ccode = [1,0,0];
    elseif ismember(c,bondinsp)
        ccode = [0,1,0];
    else
        ccode = [0,0,1];
    end
    if dx>nx/2
        dx = dx-nx;
        plot([pos(i,1),pos(i,1)+dx],[pos(i,2),pos(j,2)],'Color',ccode,'Linewidth',2);
        plot([pos(j,1)-dx,pos(j,1)],[pos(i,2),pos(j,2)],'Color',ccode,'Linewidth',2);
    elseif dx<-nx/2
        dx = dx+nx;
        plot([pos(i,1),pos(i,1)+dx],[pos(i,2),pos(j,2)],'Color',ccode,'Linewidth',2);
        plot([pos(j,1)-dx,pos(j,1)],[pos(i,2),pos(j,2)],'Color',ccode,'Linewidth',2);
    else
        plot([pos(i,1),pos(i,1)+dx],[pos(i,2),pos(j,2)],'Color',ccode,'Linewidth',2);
    end
end
plot(pos(:,1),pos(:,2),'ko','markersize',16,'markerfacecolor',ccode1,'color',ccode2,'linewidth',0.5);
xlim([-1/4,nx-1/4]);
ylim([-sqrt(3)/4,ny-sqrt(3)/4]);
%}