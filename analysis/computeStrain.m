function [strn,vor] = computeStrain(disp,pos,nb,opt)
% the function compute the strain tensor for a given displacement field
% when use the weight function by Leibler,
% opt.fun = 'weightLeibler', opt.min = 6, opt.max = 8
dim = size(pos,1);
num = size(pos,2);
eps = 0.001*max(abs(disp))/sqrt(dim);
L = max(pos(1,:))-min(pos(1,:))+0.5;

wt  = ones(size(nb,1),1); %weight function
if nargin>3 %if the weight function is determined
    wt = feval(opt.fun,pos,nb,opt);
end

if numel(disp)~=dim*num
    error('size mismatches!');
end

strn = cell(num,1);
vor  = cell(num,1);

for n = 1:num
    id1 = find(nb(:,1)==n);
    id2 = find(nb(:,2)==n);
    nj  = [nb(id1,2);nb(id2,1)]; % neighbor of n
    Amat = zeros(dim,dim);
    Dmat = zeros(dim,dim);
    for d1 = 1:dim
        dp1 = pos(d1,nj)'-pos(d1,n)';
        for ii=1:length(nj)
            if dp1(ii)>L/2
                dp1(ii) = dp1(ii)-L;
            elseif dp1(ii)<-L/2
                dp1(ii) = dp1(ii)+L;
            end
        end
        dx1 = wt([id1;id2]).*dp1;
        dxx = dx1+eps*wt([id1;id2]).*(disp(dim*(nj-1)+d1)-disp(dim*(n-1)+d1));
        for d2 = 1:dim
            dx2 = (pos(d2,nj)-pos(d2,n));
            for ii=1:length(nj)
                if dx2(ii)>L/2
                    dx2(ii) = dx2(ii)-L;
                elseif dx2(ii)<-L/2
                    dx2(ii) = dx2(ii)+L;
                end
            end
            Amat(d1,d2) = dx2*dxx;
            Dmat(d1,d2) = dx2*dx1;
        end
    end
    Fmat = Amat/Dmat;
    strn{n} = 0.5*(Fmat*Fmat'-eye(dim))/eps;
    vor{n}  = 0.5*(Fmat-Fmat')/eps;
end