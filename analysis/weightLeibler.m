function wt = weightLeibler(pos,nb,opt)
% this function computes the weight using Leibler's definition
% wt = 1 if r<opt.min
% wt = 1-(r-opt.min)/(opt.max-opt.min)
% wt = 0 if r>opt.max, not in nb
dim = size(pos,1);
%num = size(pos,2);

wt  = zeros(size(nb,1),1);
for n = 1:size(nb,1)
    n1 = nb(n,1);
    n2 = nb(n,2);
    dx = zeros(dim,1);
    for d = 1:dim
        dx(d) = pos(d,n1) - pos(d,n2); % not periodic
    end
    dr = sqrt(sum(dx.^2));
    if dr < opt.min
        wt(n) = 1;
    elseif dr<opt.max
        wt(n) = 1-(dr-opt.min)/(opt.max-opt.min);
    end
end