function pars = initialpars(num,pos,nb,ppid,L)
%ptid = find(pert~=0);
%ppid = union(ceil(ptid/2),[]);
rrid = setdiff((1:num)',ppid);
nnb  = size(nb,1);

pars.x0 = reshape(pos(rrid,:)',[],1);
pars.pos = reshape(pos(ppid,:)',[],1);
pars.L = L;

cnb=0;
cex=0;
for i=1:nnb
    n1 = nb(i,1);
    n2 = nb(i,2);
    if ismember(n1,rrid)&&ismember(n2,rrid)
        cnb=cnb+1;
        pars.nb(cnb,1) = find(rrid==n1);
        pars.nb(cnb,2) = find(rrid==n2);
        dx = pos(n1,1)-pos(n2,1);
        dy = pos(n1,2)-pos(n2,2);
        if dy>0.5*L
            dy = dy-L;
        elseif dy<-0.5*L
            dy = dy+L;
        end
        pars.nbl(cnb)  = sqrt(dx^2+dy^2);
    elseif ismember(n1,ppid)&&ismember(n2,rrid)
        cex = cex+1;
        pars.ext(cex,1) = find(ppid==n1);
        pars.ext(cex,2) = find(rrid==n2);
        dx = pos(n1,1)-pos(n2,1);
        dy = pos(n1,2)-pos(n2,2);
        if dy>0.5*L
            dy = dy-L;
        elseif dy<-0.5*L
            dy = dy+L;
        end
        pars.exl(cex)   = sqrt(dx^2+dy^2);
    elseif ismember(n1,rrid)&&ismember(n2,ppid)
        cex = cex+1;
        pars.ext(cex,1) = find(ppid==n2);
        pars.ext(cex,2) = find(rrid==n1);
        dx = pos(n1,1)-pos(n2,1);
        dy = pos(n1,2)-pos(n2,2);
        if dy>0.5*L
            dy = dy-L;
        elseif dy<-0.5*L
            dy = dy+L;
        end
        pars.exl(cex)   = sqrt(dx^2+dy^2);
    end
end