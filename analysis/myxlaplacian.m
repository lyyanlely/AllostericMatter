function [vMat,d] = myxlaplacian(Amat,q,opt)
% the function uses x-laplacian method to decompose the first q modes into
% localized ones and extended ones
if nargin<3
    opt.eta = -10.;  % negative for smallest eigenvalues
    opt.sig = 'sm';
    opt.del = 4;  % threshold for stopping remove local ones
end

dim = size(Amat,1);

%vMat = zeros(dim,q); % vectors
%d    = zeros(q,1); % eigenvalues

[V,~] = eigs(Amat,q,opt.sig);
invp  = sum(V.^4);  % inverse participation ratio
%[~,id]= max(invp);
%vMat(:,1) = V(:,id);
%d(1) = D(id,id);

Lmat = Amat-opt.eta*diag(sum(V(:,invp>opt.del).^2,2));

[vMat,D] = eigs(Lmat,q,opt.sig);

d = diag(D);