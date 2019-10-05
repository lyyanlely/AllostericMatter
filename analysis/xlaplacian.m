function vMat = xlaplacian(Amat,q,opt)
% the function uses x-laplacian method to decompose the first q modes into
% localized ones and extended ones
if nargin<3
    opt.ntr = 6;  % number of trivial modes
    opt.eta = -.2;  % negative for smallest eigenvalues
    opt.sig = 'sm';
    opt.del = 4;  % threshold for stopping remove local ones
end

dim = size(Amat,1);

vMat = zeros(dim,q); % vectors
%d    = zeros(q,1); % eigenvalues

[V,D] = eigs(Amat,q+opt.ntr,opt.sig);
id  = abs(diag(D))>1e-8; %
Vtr = V(:,~id); % trivial modes
Pnr = eye(dim)-Vtr*Vtr';  % projection operator
V   = V(:,id); % nontrivial modes
invp  = sum(V.^4);  % inverse participation ratio
[~,id]= max(invp);
vMat(:,1) = V(:,id);
%d(1) = D(id,id);

Lmat = Amat-opt.eta*diag(V(:,id));

for i = 2:q
    LPmat = Pnr*Lmat*Pnr;
    [V,D] = eigs(LPmat,q-i+1+opt.ntr,opt.sig);
    id  = abs(diag(D))>1e-8; %
    V   = V(:,id); % nontrivial modes
    invp  = sum(V.^4);  % inverse participation ratio
    if any(invp*dim>opt.del)
        [~,id]= max(invp);
        vMat(:,i) = V(:,id);
        %d(i) = D(id,id);
        Lmat = Lmat-opt.eta*diag(V(:,id).^2);
    else
        vMat(:,i:q) = V;
        %d(i:q) = diag(D);
        break;
    end
end