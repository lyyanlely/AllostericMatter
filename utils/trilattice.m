function [pos,bondnb,Smat,bpos]=trilattice(xl, yl, eps)
% build a network of xl and yl size eps distortion
% xl    = 20;
% yl    = 20;
% eps   = 0.2;   % distort size
% pos are position of particles
% bondnb gives the indices of the two particles a bond connecting in each
% row
% Smat gives the structural matrix
% bpos gives the midpoint position of a bond
warning('off','all');
%% size of system and size dependent parameters
num   = xl*yl;
dim   = 2;
xbound = xl;
ybound = yl*sqrt(3)/2;
%nmax  = 3*num-2*xl;

%% construct the embedded network
bondnb = zeros(3*num,2);
nnb    = 0;
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
            posx(n) = posx(n);% + eps*sqrt(3)/2;
            posy(n) = posy(n);% + eps/2;
        end
    else
        if mod(nx+mod(ny/2-1,2),2)
            posx(n) = posx(n);
            posy(n) = posy(n);% - eps;
        else
            posx(n) = posx(n);% - eps*sqrt(3)/2;
            posy(n) = posy(n);%% + eps/2;
        end
    end
    %}
    theta = 2*pi*rand(1);
    posx(n) = posx(n) + eps*cos(theta);
    posy(n) = posy(n) + eps*sin(theta);
    %if nx<xl
    nnb = nnb+1;
    bondnb(nnb,1) = n;
    bondnb(nnb,2) = n+1-xl*(nx==xl);  % periodic in horizontal direction
    %end
    if ny<yl  % the top line is not periodic
        %if nx<xl||mod(ny,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl+(1-xl*(nx==xl))*mod(ny-1,2)-num*(ny==yl); %
        %end
        %if nx>1||mod(ny+1,2)
        nnb = nnb+1;
        bondnb(nnb,1) = n;
        bondnb(nnb,2) = n+xl-(1-xl*(nx==1))*mod(ny,2)-num*(ny==yl);  %
        %end
    end
end
bondnb = bondnb(1:nnb,:);
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
    %
    if dely(n)>ybound/2  % not periodic
        dely(n) = dely(n) - ybound;
    elseif dely(n)<-ybound/2
        dely(n) = dely(n) + ybound;
    end
    %}
    bpx(n) = posx(jj)+delx(n)/2;
    bpy(n) = posy(jj)+dely(n)/2;
    delr(n) = sqrt(delx(n)^2+dely(n)^2);
    idx(4*(n-1)+1) = n;
    jdx(4*(n-1)+1) = 2*(ii-1)+1;    % 2
    val(4*(n-1)+1) = delx(n)/delr(n);
    idx(4*(n-1)+2) = n;
    jdx(4*(n-1)+2) = 2*(ii-1)+2;    % 1
    val(4*(n-1)+2) = dely(n)/delr(n);
    idx(4*(n-1)+3) = n;
    jdx(4*(n-1)+3) = 2*(jj-1)+1;    % 2
    val(4*(n-1)+3) = -delx(n)/delr(n);
    idx(4*(n-1)+4) = n;
    jdx(4*(n-1)+4) = 2*(jj-1)+2;    % 1
    val(4*(n-1)+4) = -dely(n)/delr(n);
end
Smat = sparse(idx,jdx,val,nnb,dim*num);
pos  = [posx/xbound;posy/xbound;]';  % x y
bpos = [bpx/xbound;bpy/xbound;]';    % x y