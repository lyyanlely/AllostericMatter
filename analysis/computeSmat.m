function Smat=computeSmat(pos,nb,pflg)
nnb  = size(nb,1);
num  = size(pos,1);
dim  = size(pos,2);
xbound = ceil(max(pos(:,1)));

delx = zeros(1,nnb);
dely = zeros(1,nnb);
delr = zeros(1,nnb);
idx = zeros(1,4*nnb);
jdx = zeros(1,4*nnb);
val = zeros(1,4*nnb);
for n = 1:nnb
    ii = nb(n,1);
    jj = nb(n,2);
    delx(n) = pos(ii,1) - pos(jj,1);
    dely(n) = pos(ii,2) - pos(jj,2);
    if pflg
        if delx(n)>xbound/2   % periodic
            delx(n) = delx(n) - xbound;
        elseif delx(n)<-xbound/2
            delx(n) = delx(n) + xbound;
        end
    end
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
Smat = sparse(idx,jdx,val,nnb,dim*num);