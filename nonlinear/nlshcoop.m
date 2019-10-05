function [ecoop,er,eb,eb2] = nlshcoop(kk,a,b0,k1,opts)
% e_w=0.5*(0.5*dl^2-a^2)*dl^2, dl=l/l0-1 normalized by the lattice constant
if nargin<5
    opts = [];
end
if isfield(opts,'gifname')
    gifname = opts.gifname;
else
    gifname = 'nlshear_bind_b3.gif';
end
if isfield(opts,'num')
    num = opts.num;
else
    num = 400;   % size of the system
end
if isfield(opts,'prtlevel')
    options.prtlevel = opts.prtlevel;
else
    options.prtlevel = 0;
end
L = sqrt(num);  
dim = 2;   % dimension
reps=0.01/L;   % random displacement before conjugate gradient

upp = 0.5+0.05/L; % right boundary
low = 0.5-0.75/L; % left boundary in between are weak springs
x0 = 0.5-.4/L;

%k1 = 10.;
%kk = .05; % stiffness of the weak springs
%a  = .5;  % mismatch between the lattice constant l0 and the restlengths of weak springs
tt = -((2*a+0.5)/L):0.1/L:((2*a+0.5)/L); % angle along the shear mode
lb = length(tt);
%b0 = 0.;  % reference extenal field on weak springs
%tid  = 10;

% build the network on triangular lattice
[pos,nb,~,bpos]  = trilattice(sqrt(num),sqrt(num),0.);
% make the bonds in [low,upp] softer
flg = bpos(:,1)>low&bpos(:,1)<upp; 
bids = nb(find(flg,1),:);  % id of binding sites
vbid = [2*bids(1)-1,2*bids(1),2*bids(2)-1,2*bids(2)];
%nbid = setdiff(1:dim*num,vbid);
bidt = nb(find(flg,1,'last'),:);  % id of binding sites
vbit = [2*bidt(1)-1,2*bidt(1),2*bidt(2)-1,2*bidt(2)];
%nbit = setdiff(1:dim*num,vbit);

pos0 = pos;
ts = zeros(dim*num,2);
ts(1:2:end,1) = 1;
ts(:,1) = ts(:,1)/norm(ts(:,1));
ts(2:2:end,2) = 1;
ts(:,2) = ts(:,2)/norm(ts(:,2));


e0 = zeros(1,lb);     % energy before relaxation
for i=1:lb
    pars = initialparsBE(pos0,nb,flg,L,kk,a,0); % initialize parameters for CG
    pars.nvar  = dim*num;       % number of variables
    pars.fgname= 'BistableE'; %'BisEnRot'; %  % function of computing energy to be minimized
    [x,~] = shearsquare(pars.x0,tt(i),L,x0); % rotate the block to the desired angle along the shear mode
    [e,~] = BistableE(x,pars);  % compute the corresponding elastic energy before relaxation
    e0(i) = e;
    if i==10
        x2 = x;
    elseif i==22
        x1 = x;
    end
end
[~,lcs] = findpeaks(-e0);
tid = lcs(1);
rid = lcs(2);
%er = e0(rid);  % reference energy

t0 = tt(tid);    % binding theta
%
pars = initialparsBE(pos0,nb,flg,L,kk,a,b0); % initialize parameters for CG
pars.nvar  = dim*num;       % number of variables
pars.fgname= 'BistableE'; %'BisEnRot'; %  % function of computing energy to be minimized
[x,v] = shearsquare(pars.x0,t0,L,x0); % rotate the block to the desired angle along the shear mode
v = v-ts(:,1)'*v*ts(:,1)-ts(:,2)'*v*ts(:,2);
v = v/norm(v);
options.x0 = x+reps*(rand(size(pars.x0))-0.5);  % add some random perturbation to the current location to break symmetry in triangular lattice
options.maxit = 10000;
pars.v = v; 
pars.mid = x0;
[x,e,~,~,~]=nlcg(pars,options);
pars = rmfield(pars,'v');
er = e;
%bp0 = x(vbid); % - pars.x0(vbid);
% constraint at the stimulus
n1 = bids(1);
n2 = bids(2);
dx = x(2*n2-1)-x(2*n1-1);
dy = x(2*n2)-x(2*n1);
dr0 = sqrt(dx^2+dy^2);
th0 = atan(dy/dx);
if dx<0
    th0 = th0+pi;
end
% constraint at the target
n1 = bidt(1);
n2 = bidt(2);
dx = x(2*n2-1)-x(2*n1-1);
dy = x(2*n2)-x(2*n1);
dr1 = sqrt(dx^2+dy^2);
th1 = atan(dy/dx);
if dx<0
    th1 = th1+pi;
end
%}
%ee = zeros(1,lb);     % energy after relaxation (no relax along shear mode)
%eb = zeros(1,lb);     % energy with binding
%eb2 = zeros(1,lb);     % energy with binding both

pars.mid = x0;
options.x0 = x1+reps*(rand(size(pars.x0))-0.5);  % add some random perturbation to the current location to break symmetry in triangular lattice
[x,e,~,~,~]=nlcg(pars,options);   % use the nonlinear conjugate gradient method to relax
ee = e;
options.x0 = x2+reps*(rand(size(pars.x0))-0.5);  % add some random perturbation to the current location to break symmetry in triangular lattice
[x,e,~,~,~]=nlcg(pars,options);   % use the nonlinear conjugate gradient method to relax
if ee>e
    ee = e;
elseif abs(ee-e)<1e-10
    ecoop = 0;
    eb = e;
    eb2 = e;
    return
end
% energy with binding
pars.bid = bids;  % locations the motion is zero
pars.k0  = k1;
pars.th0 = th0;
pars.r0  = dr0;
options.x0 = x1+reps*(rand(size(pars.x0))-0.5);
%options.x0(vbid) = bp0;
[x,e,~,~,~]=nlcg(pars,options);
n1 = bids(1);
n2 = bids(2);
dx = x(2*n2-1)-x(2*n1-1);
dy = x(2*n2)-x(2*n1);
dr = sqrt(dx^2+dy^2);
th = atan(dy/dx);
if dx<0
    th = th+pi;
end
eb = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2;
options.x0 = x2+reps*(rand(size(pars.x0))-0.5);
%options.x0(vbid) = bp0;
[x,e,~,~,~]=nlcg(pars,options);
if eb>e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2
    eb = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2;
end

% energy with binding both
pars.bit = bidt;  % locations the motion is zero
pars.k0  = k1;
pars.th1 = th1;
pars.r1  = dr1;
options.x0 = x1+reps*(rand(size(pars.x0))-0.5);
[x,e,~,~,~]=nlcg(pars,options);
n1 = bids(1);
n2 = bids(2);
dx = x(2*n2-1)-x(2*n1-1);
dy = x(2*n2)-x(2*n1);
dr = sqrt(dx^2+dy^2);
th = atan(dy/dx);
if dx<0
    th = th+pi;
end
n1 = bidt(1);
n2 = bidt(2);
dx = x(2*n2-1)-x(2*n1-1);
dy = x(2*n2)-x(2*n1);
dr = sqrt(dx^2+dy^2);
th = atan(dy/dx);
if dx<0
    th = th+pi;
end
eb2 = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2-0.5*k1*(dr-dr1)^2-0.5*k1*(th-th1)^2*pars.l0^2;
options.x0 = x2+reps*(rand(size(pars.x0))-0.5);
[x,e,~,~,~]=nlcg(pars,options);
if eb2>e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2-0.5*k1*(dr-dr1)^2-0.5*k1*(th-th1)^2*pars.l0^2
    eb2 = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2-0.5*k1*(dr-dr1)^2-0.5*k1*(th-th1)^2*pars.l0^2;
end

ecoop = 2*eb - eb2 - ee;