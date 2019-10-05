%function plotlever()
% this function generates a figure for lever
%num   = 256;
lx    = 16;
ly    = 16;
dim   = 2;
num   = lx*ly;
z     = 4;
x0    = 0.5;
r     = 4;
vi    = 1;
flr   = 0; % 0 exct l resp r 1 exct r resp l
%% color
respc = [0.,0.6,0.8];
exctc = [0.7,0.2,0.7];
netwc = [0.1,0.1,0.1];
%% network
[pos,nb,fl,fr]  = cutNetwork(lx,ly,dim,z,r,x0);
[dl,dr,Tlr,Trl] = LRCoupling(lx,ly,dim,z,x0,r);
[Smat,~] = StructureMat(pos,nb);
pos(:,1) = pos(:,1)-x0;
pos(pos(:,1)<0,1) = pos(pos(:,1)<0,1)+lx/ly;
%% magnifi mode
if flr>0
    [U,S,V] = svd(Tlr);  % singular value decomposition U*S=Tlr*V
    ie = fr;
    ir = fl;
else
    [U,S,V] = svd(Trl);
    ie = fl;
    ir = fr;
end
%% find coupling field
Z  = null(full(Smat'*Smat));
v1 = zeros(dim*num,1);
v2 = zeros(dim*num,1);
v1(1:dim:end) = 1/sqrt(num);
v2(2:dim:end) = 1/sqrt(num);
Zt = Z-v1*(v1'*Z)-v2*(v2'*Z);
Zn = orth(Zt);
Zn = Zn(:,1:end-2);
Il = zeros(length(ie)*2,numel(pos));
Ir = zeros(length(ir)*2,numel(pos));
for i=1:length(ie)
    Il(2*i-1,2*ie(i)-1) = 1;
    Il(2*i,2*ie(i)) = 1;
end
for i=1:length(ir)
    Ir(2*i-1,2*ir(i)-1) = 1;
    Ir(2*i,2*ir(i)) = 1;
end
plr = [Il;Ir];
pZ  = plr*Zn;

vecb = [V(:,vi);S(vi,vi)*U(:,vi)];
coef = pZ\vecb;
fmod = Zn*coef; %/S(vi,vi)*sc;

plotNetwork0(pos,nb,lx/ly)
quiver(pos(ie,1),pos(ie,2),V(1:2:end,vi),V(2:2:end,vi),'linewidth',2,'AutoScale','on','color',exctc)
quiver(pos(ir,1),pos(ir,2),U(1:2:end,vi),U(2:2:end,vi),'linewidth',2,'AutoScale','on','color',respc)
quiver(pos(:,1),pos(:,2),fmod(1:2:end),fmod(2:2:end),'linewidth',1,'AutoScale','on','color',netwc)