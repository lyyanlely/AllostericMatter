function [strn,vor,bk,sh,vr] = computeStrain_real(pos1,pos2,nb,opt)
% the function compute the strain tensor for a given displacement field
% when use the weight function by Leibler,
% opt.fun = 'weightLeibler', opt.min = 6, opt.max = 8
dim = size(pos1,1);
num = size(pos1,2);
sample=20;
%eps = 0.001*max(abs(disp))/sqrt(dim);
eps=1;
wt  = ones(size(nb,1),1); %weight function
if nargin>3 %if the weight function is determined
    wt = feval(opt.fun,pos1,nb,opt);
end

strn = cell(num,1);
vor  = cell(num,1);

bk  = zeros(num,1);
sh  = zeros(num,1);
vr  = zeros(num,1);

for i = 1:sample
for n = 1:num
    id1 = find(nb(:,1)==n);
    id2 = find(nb(:,2)==n); % neighbor of n
    nj  = [nb(id1,2);nb(id2,1)];
    Amat = zeros(dim,dim);
    Dmat = zeros(dim,dim);
    for d1 = 1:dim
        dx1 = wt([id1;id2]).*(pos1(d1,nj)'-pos1(d1,n)');
        dxx = eps*wt([id1;id2]).*(pos2(d1,num*(i-1)+nj)'-pos2(d1,num*(i-1)+n)');
        for d2 = 1:dim
            dx2 = (pos1(d2,nj)-pos1(d2,n));
            Amat(d1,d2) = dx2*dxx;
            Dmat(d1,d2) = dx2*dx1;
        end
    end
    
    Fmat = Amat/Dmat;
    strn{n} = 0.5*(Fmat*Fmat'-eye(dim))/eps;
    vor{n}  = 0.5*(Fmat-Fmat')/eps;
    sttr = trace(strn{n});
    bk(n) = bk(n) + sum(diag(strn{n}).^2);
    sh(n) = sh(n) + 0.5*sum(sum((strn{n}-sttr*eye(dim)).^2));
    vr(n) = vr(n) + vor{n}(1,2);    
    %bk(n) = sum(diag(strn{n}).^2);
    %sh(n) = 0.5*sum(sum((strn{n}-sttr*eye(dim)).^2));
    %vr(n) = vor{n}(1,2);
end
end
bk=bk/sample;
sh=sh/sample;
vr=vr/sample;