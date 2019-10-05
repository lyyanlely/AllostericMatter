function vMat = xlaplacian1(Amat,q,opt)
% the function uses x-laplacian method to decompose the first q modes into
% localized ones and extended ones
if nargin<3
    opt.dim = 3;  % spatial dimension
    opt.ntr = 6;  % number of trivial modes
    opt.eta = -1.;  % negative for smallest eigenvalues
    opt.sig = 'sm';
    opt.del = 12;  % threshold for stopping remove local ones
end
%q = 1;

dim = size(Amat,1);

%vMat = zeros(dim,q); % vectors
%d    = zeros(q,1); % eigenvalues

[V,D] = eigs(Amat,q+opt.ntr,opt.sig);
id  = abs(diag(D))>1e-8; %
Vtr = V(:,~id); % trivial modes
Pnr = eye(dim)-Vtr*Vtr';  % projection operator
V   = V(:,id); % nontrivial modes
invp  = sum(V.^4);  % inverse participation ratio
%[~,id]= max(invp);
%vMat = V(:,id)';
vMat = [];
%d(1) = D(id,id);
Lmat = Amat;

while any(dim*invp>opt.del)
    [~,id]= max(invp);
    if isempty(vMat)
        vMat = V(:,id)';
    else
        vMat  = [vMat;V(:,id)'];
    end
    Lmat = Lmat-opt.eta*diag(V(:,id).^2);
    LPmat = Pnr*Lmat*Pnr;
    [V,D] = eigs(LPmat,q+opt.ntr,opt.sig);
    id  = abs(diag(D))>1e-8; %
    V   = V(:,id); % nontrivial modes
    invp  = sum(V.^4);  % inverse participation ratio
end

vMat = [vMat;V'];
vMat = vMat';