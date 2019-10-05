function sites = findIsostatic(sites,boundary,bondnb,pos)
% this function identifies the maximum isostatic network surrounding the
% initial sites (initially rigid networks)
dim = 2;
f   = .5; % fraction number counted for the external sites
num = max(max(bondnb)); % total number of sites
neighb = cell(num,1);   % neighbor sites for each site
nnb = zeros(num,1);     % number of neighbors
dist = zeros(num,1);    % distance to the initial sites

p0 = [mean(pos(sites,1)),mean(pos(sites,2))];
for i = 1:num
    flg = sum(ismember(bondnb,i),2)>0;
    neighb{i} = setdiff(reshape(bondnb(flg,:),1,[]),i);
    nnb(i)    = length(neighb{i});
    dx = pos(i,1)-p0(1);
    if dx>0.5*round(max(pos(:,1)))
        dx = dx-round(max(pos(:,1)));
    elseif dx<-0.5*round(max(pos(:,1)))
        dx = dx+round(max(pos(:,1)));
    end
    dy = pos(i,2)-p0(2);
    dist(i) = sqrt(dx^2+dy^2);
end

%intb  = []; % internal bonds
%extb  = []; % external bonds
nset  = sum(ismember(bondnb,sites),2);
intb  = find(nset==2);
extb  = find(nset==1);

nsite = length(sites);  % number of sites
nbond = length(intb)+f*length(extb); %0.5*sum(nnb(sites));  % number of shared bonds
if nbond>nsite*dim   %<
    sites = [];
    %print('Initial sites are not in an isostatic cluster!');
    return
else
    u = setdiff(1:num,boundary'); % set of sites not tested
    t = [];  % set of sites to be tested
    for i = 1:length(sites)
        t = [t,neighb{sites(i)}];
    end
    t = intersect(t,u);
    [n0,id] = min(nnb(t));  %max
    if sum(nnb(t)==n0)>1
        flg = find(nnb(t)==n0);
        [~,i1] = min(dist(t(flg)));
        id  = flg(i1);
    end
    t1 = t(id);
    t(id) = [];
    u  = setdiff(u,t1);
    nset  = sum(ismember(bondnb,[sites,t1]),2);
    intb  = find(nset==2);
    extb  = find(nset==1);
    while length(intb)+f*length(extb)<=(nsite+1)*dim % floppy not makes the network rigid >
        sites = [sites,t1];
        nsite = nsite+1;
        %nbond = length(ibtmp)+f*length(ebtmp);
        t = union(t,intersect(neighb{t1},u)); % to test
        tmp = t;
        while ~isempty(tmp)
            t1 = tmp(1);
            tmp(1) = [];
            if all(ismember(neighb{t1},sites))
                t(t==t1) = [];
                u     = setdiff(u,t1);
                sites = [sites,t1];
                nsite = nsite+1;
            end
        end
        [n0,id] = min(nnb(t));  %max
        if sum(nnb(t)==n0)>1
            flg = find(nnb(t)==n0);
            [~,i1] = min(dist(t(flg)));
            id  = flg(i1);
        end
        t1 = t(id);
        t(id) = [];
        u  = setdiff(u,t1);
        nset  = sum(ismember(bondnb,[sites,t1]),2);
        intb  = find(nset==2);
        extb  = find(nset==1);
    end
end