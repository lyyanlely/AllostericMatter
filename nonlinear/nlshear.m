function [ecoop,ee,eb,eb2] = nlshear(kk,a,b0,k1,opts)
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
    num = 100;   % size of the system
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
end
[~,lcs] = findpeaks(-e0);
tid = lcs(1);
rid = lcs(2);

t0 = tt(tid);    % binding theta
pars = initialparsBE(pos0,nb,flg,L,kk,a,b0); % initialize parameters for CG
pars.nvar  = dim*num;       % number of variables
pars.fgname= 'BistableE'; %'BisEnRot'; %  % function of computing energy to be minimized
%{
[x,v] = shearsquare(pars.x0,t0,L,x0); % rotate the block to the desired angle along the shear mode
v = v-ts(:,1)'*v*ts(:,1)-ts(:,2)'*v*ts(:,2);
v = v/norm(v);
options.x0 = x+reps*(rand(size(pars.x0))-0.5);  % add some random perturbation to the current location to break symmetry in triangular lattice
options.maxit = 10000;
pars.v = v; 
pars.mid = x0;
[x,e,~,~,~]=nlcg(pars,options);
pars = rmfield(pars,'v');
%}
% get the rest lengths of two global states
nnb  = size(nb,1);
pars.x1 = zeros(nnb,2);  % two restlengths that are stable
[x,~] = shearsquare(pars.x0,-t0,L,x0);
for fg = [1,-1,2]
    id = find(pars.flag==fg,1);
    n1 = nb(id,1);
    n2 = nb(id,2);
    dr = sqrt((x(2*n2-1)-x(2*n1-1))^2+(x(2*n2)-x(2*n1))^2);
    pars.x1(pars.flag==fg,1) = dr;
end
[x,~] = shearsquare(pars.x0,t0,L,x0);
for fg = [1,-1,2]
    id = find(pars.flag==fg,1);
    n1 = nb(id,1);
    n2 = nb(id,2);
    dr = sqrt((x(2*n2-1)-x(2*n1-1))^2+(x(2*n2)-x(2*n1))^2);
    pars.x1(pars.flag==fg,2) = dr;
end
l0 = pars.l0/pars.a;
pars.p = zeros(nnb,5);
for fg = [1,-1]
    id = find(pars.flag==fg,1);
    matA = zeros(5);
    vecB = zeros(5,1);
    matA(1,1) = 1; 
    vecB(1)   = kk/4*pars.a^2;
    matA(2,:) = [l0^4,l0^3,l0^2,l0,1];
    matA(3,:) = [4*pars.x1(id,1)^3/pars.a^3,3*pars.x1(id,1)^2/pars.a^2,2*pars.x1(id,1)/pars.a,1,0];
    matA(4,:) = [4*pars.x1(id,2)^3/pars.a^3,3*pars.x1(id,2)^2/pars.a^2,2*pars.x1(id,2)/pars.a,1,0];
    matA(5,:) = [(pars.x1(id,1)^4-pars.x1(id,2)^4)/pars.a^4,(pars.x1(id,1)^3-pars.x1(id,2)^3)/pars.a^3,(pars.x1(id,1)^2-pars.x1(id,2)^2)/pars.a^2,(pars.x1(id,1)-pars.x1(id,2))/pars.a,0];
    vecB(5)   = -kk*b0*pars.a^2;
    pars.p(pars.flag==fg,:) = repmat((matA\vecB)',sum(pars.flag==fg),1);
end
fg = 2;
id = find(pars.flag==fg,1);
a0   = kk/2*pars.a^2;
pars.p(pars.flag==fg,3:5) = repmat([a0,-2*a0*pars.x1(id,1)/pars.a,a0*(pars.x1(id,1)/pars.a)^2],sum(pars.flag==fg),1);
para = pars.p;
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

ee = zeros(1,lb);     % energy after relaxation (no relax along shear mode)
eb = zeros(1,lb);     % energy with binding
eb2 = zeros(1,lb);     % energy with binding both

for i=1:lb
    if i==ceil(lb/2)
        1;
    end
    pars = initialparsBE(pos0,nb,flg,L,kk,a,b0); % initialize parameters for CG
    pars.p = para;
    pars.nvar  = dim*num;       % number of variables
    pars.fgname= 'BistableE'; %'BisEnRot'; %  % function of computing energy to be minimized
    [x,v] = shearsquare(pars.x0,tt(i),L,x0); % rotate the block to the desired angle along the shear mode
    v = v-ts(:,1)'*v*ts(:,1)-ts(:,2)'*v*ts(:,2);
    v = v/norm(v);
    %[x,~,~,~,~]=nlcg(pars,options);
    %mm = BistableMat(x,pars);   % compute the dynamical matrix for current position and bistable network energy function
    %[V,d]=eigs(mm,3,'sm');      % find the three modes of smallest absolute value of eigenenergy (two are translations)
    %id = find(abs(diag(d))>1e-8);    % the other is the shear mode (mechanism) we don't want to relax along
    %pars.b = b0;
    %[e,~] = BistableE(x,pars);  % compute the corresponding elastic energy before relaxation
    %e0(i) = e;
    options.x0 = x+reps*(rand(size(pars.x0))-0.5);  % add some random perturbation to the current location to break symmetry in triangular lattice
    options.maxit = 10000;
    pars.v = v; %V(:,id); % the shear mode (mechanism) we don't relax along
    pars.mid = x0;
    [x,e,~,~,~]=nlcg(pars,options);   % use the nonlinear conjugate gradient method to relax
    ee(i) = e;
    % energy with binding
    pars.bid = bids;  % locations the motion is zero
    pars.k0  = k1;
    pars.th0 = th0;
    pars.r0  = dr0;
    options.x0 = x+reps*(rand(size(pars.x0))-0.5);
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
    eb(i) = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2;
    
    % energy with binding both
    pars.bit = bidt;  % locations the motion is zero
    pars.k0  = k1;
    pars.th1 = th1;
    pars.r1  = dr1;
    options.x0 = x+reps*(rand(size(pars.x0))-0.5);
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
    eb2(i) = e-0.5*k1*(dr-dr0)^2-0.5*k1*(th-th0)^2*pars.l0^2;
    n1 = bidt(1);
    n2 = bidt(2);
    dx = x(2*n2-1)-x(2*n1-1);
    dy = x(2*n2)-x(2*n1);
    dr = sqrt(dx^2+dy^2);
    th = atan(dy/dx);
    if dx<0
        th = th+pi;
    end
    eb2(i) = eb2(i)-0.5*k1*(dr-dr1)^2-0.5*k1*(th-th1)^2*pars.l0^2;
    pars = rmfield(pars,{'bid','bit'});
    pars = rmfield(pars,'v');   % remove the shear mode from the parameters for the next round of computation
    % plot the relaxed structures and save them in a gif file
    %{
    pos = reshape(x,2,[])';
    plotNetwork0(pos,nb,pars.xb,flg);
    f  = getframe;
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    if i==1
        imwrite(im,map,gifname,'gif','DelayTime',0,'LoopCount',inf);
    else
        imwrite(im,map,gifname,'gif','WriteMode','append');
    end
    close all
    %}
end

%[e2,id] = min(eb2);
%em = min(ee);
%ecoop = 2*min(eb) - e2 - em;
%
[ems,lcs] = findpeaks(-ee);
if length(lcs)<2 || ~ismember(lcs(1),tid-1:tid+1)
    ecoop = 0;
else
    [e2,id] = min(eb2);
    em = min(ee); %it = lcs(2);
    if ismember(id,tid-1:tid+1) %&& ismember(it,rid-1:rid+1)
        ecoop = 2*min(eb) - e2 - em;
    else
        ecoop = 0;
    end
end
%}