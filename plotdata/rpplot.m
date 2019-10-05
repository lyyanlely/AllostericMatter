function rpplot(replica,config,num,linnum,numcon,bondnb)
hlayer = zeros(linnum,2*linnum);
vlayer= zeros(linnum,2*linnum);
for ii=1:linnum
    for jj=1:linnum
        hlayer(ii, 2*(jj-1)+1) = 3*linnum*(ii-1)+3*(jj-1)+2;
        hlayer(ii, 2*jj)       = 3*linnum*(ii-1)+3*jj;
        vlayer(ii, 2*(jj-1)+1) = 3*linnum*(jj-1)+3*(ii-1)+1;
        vlayer(ii, 2*jj)       = 3*linnum*(jj-1)+3*(ii-1)+2+4*mod(jj,2)-mod(jj,2)*floor(ii/linnum)*3*linnum;
    end
end
state  = config(replica,:);
occupy = find(state>0);
pebble = zeros(2,num);
block  = -2*ones(1,num);
btree  = zeros(1,num);
shell  = zeros(1,num);
jmax   = 0;
ndep   = 0;
tag    = 0;
stresd = zeros(1,numcon);
coverindex = zeros(1,numcon);  %a number labels a pebble from which site covers, 0 no covering pebble overconstraint bond
% adding bonds through pebble game
for n=1:numcon
    newbd = bondnb(occupy(n),:);
    [ndep_temp,tag,pebble,shell,block,btree,jmax] = check(newbd(1),newbd(2),ndep,tag,pebble,shell,block,btree,jmax);
    if ndep_temp>ndep
        % when the bond is redundant, all bonds along the path of finding
        % free pebble is in one overconstraint region.
        pbpath = shell(1:jmax);
        ocindx = setdiff(stresd,0);
        if isempty(ocindx)==0
            for tags=ocindx
                mbds  = find(stresd==tags); %marked bonds
                msts  = union(reshape(bondnb(occupy(mbds),:),1,2*length(mbds)),[]); %marked sites
                if length(intersect(msts,pbpath))>2
                    pbpath = union(pbpath,msts);
                end
            end
        end
        bonds  = ismember(occupy,find(ismember(ismember(bondnb,pbpath),[1 1],'rows')));
        stresd(bonds)=tag;
     
        coverindex(n) = 1;
        ndep = ndep_temp;
    end
end
mark     = zeros(1,numcon);
nrigid   = zeros(1,numcon);
rigidset = zeros(numcon);
[~,nrigid,rigidset] = idRigid(mark,nrigid,rigidset,numcon,bondnb(occupy,:),tag,pebble,shell,block,btree);

% find pivot sites
        % spanning: a path appears in every layer horizontally or vertically.
spanh = zeros(1,numcon);
spanv = zeros(1,numcon);
span  = zeros(1,numcon);
for n = 1:numcon
    rclust = occupy(rigidset(n,1:nrigid(n)));
    if nrigid(n)>=linnum
        spanh(n) = 1-isempty(intersect(rclust,hlayer(1,:)));
        spanv(n) = 1-isempty(intersect(rclust,vlayer(1,:)));
        for i=2:linnum
            spanh(n) = spanh(n)*(1-isempty(intersect(rclust,hlayer(i,:))));
            spanv(n) = spanv(n)*(1-isempty(intersect(rclust,vlayer(i,:))));
        end
        span(n)  = spanh(n)+spanv(n)-spanh(n)*spanv(n);
    end
end
spanclust= span==1;
stressb=find(stresd>0);
sum(state(occupy(stressb))>1)
sum(state>1)
sum(state(occupy(stressb))<1)
sum(state<1&state>0)
x = zeros(1,num);
y = zeros(1,num);
for n=1:num
    ny = ceil(n/linnum);
    nx = n - (ny-1)*linnum;
    x(n) = nx-1+mod(ny-1,2)/2;
    y(n) = sqrt(3)/2*ny - sqrt(3)/4;
end
for n=1:numcon
    ns=occupy(n);
    ip=ceil(ns/3);
    i =ns-3*(ip-1);
    %{
    if ismember(n,stressb)
        colorc=[0,54,45]/255;
    elseif ismember(n,rigidset(spanclust,:))
        colorc=[101,248,177]/255;
    else
        colorc=[5,148,229]/255;
    end
    %}
    colorc = state(ns)/2;
    if i==1
        dx=1;
        dy=0;
    elseif i==2
        dx=0.5;
        dy=sqrt(3)/2;
    else
        dx=-0.5;
        dy=sqrt(3)/2;
    end
    patch([x(ip),x(ip)+dx],[y(ip),y(ip)+dy],colorc*ones(1,2),'Edgecolor','interp','linewidth',2);
    hold all
end
axis equal off
xlim([1.5,14.5]);
ylim([1,14]);
set(gca,'xtick',[]);set(gca,'ytick',[]); box on;