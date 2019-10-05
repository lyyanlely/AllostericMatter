function [Smatrix,bpos] = SMfixV(pos,bondnb,nfix)
% this function compute the structure matrix
nnb  = size(bondnb,1);
num  = size(pos,1);
dim  = size(pos,2);
nf   = length(nfix);
mob  = setdiff(1:num,nfix); % set of mobile particles
nm   = num-nf;
posx = pos(:,1);
posy = pos(:,2);

xbound = mp(1);
ybound = mp(1);

bpx  = mp(zeros(1,nnb));
bpy  = mp(zeros(1,nnb));
delx = mp(zeros(1,nnb));
dely = mp(zeros(1,nnb));
delr = mp(zeros(1,nnb));
idx = zeros(1,4*nnb);
jdx = zeros(1,4*nnb);
val = mp(zeros(1,4*nnb));
for n = 1:nnb
    ii = bondnb(n,1);
    jj = bondnb(n,2);
    delx(n) = posx(ii) - posx(jj);
    dely(n) = posy(ii) - posy(jj);
    i0 = find(ismember(mob,ii));
    j0 = find(ismember(mob,jj));
    while delx(n)>xbound/2   % periodic
        delx(n) = delx(n) - xbound;
    end
    while delx(n)<-xbound/2
        delx(n) = delx(n) + xbound;
    end
    bpx(n) = posx(ii)+delx(n)/2;
    bpy(n) = posy(ii)+dely(n)/2;
    if dely(n)>ybound/2  % not periodic
        dely(n) = dely(n) - ybound;
    elseif dely(n)<-ybound/2
        dely(n) = dely(n) + ybound;
    end
    delr(n) = sqrt(delx(n)^2+dely(n)^2);
    if ~isempty(i0)
        idx(4*(n-1)+1) = n;
        jdx(4*(n-1)+1) = 2*(i0-1)+1;
        val(4*(n-1)+1) = delx(n)/delr(n);
        idx(4*(n-1)+2) = n;
        jdx(4*(n-1)+2) = 2*(i0-1)+2;
        val(4*(n-1)+2) = dely(n)/delr(n);
    end
    if ~isempty(j0)
        idx(4*(n-1)+3) = n;
        jdx(4*(n-1)+3) = 2*(j0-1)+1;
        val(4*(n-1)+3) = -delx(n)/delr(n);
        idx(4*(n-1)+4) = n;
        jdx(4*(n-1)+4) = 2*(j0-1)+2;
        val(4*(n-1)+4) = -dely(n)/delr(n);
    end
end
Smatrix = mp(zeros(nnb,dim*nm));
for n = 1:4*nnb
    if idx(n)>0
        Smatrix(idx(n),jdx(n)) = val(n);
    end
end
%Smatrix = sparse(idx(idx>0),jdx(idx>0),val(idx>0),nnb,dim*nm);
bpos    = [bpx;bpy]';