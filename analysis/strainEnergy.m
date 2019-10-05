function [bk,sh] = strainEnergy(strn)
% this function compute the local bulk and shear energy from the strain
% tensor
num = size(strn,1);
dim = size(strn{1},1);

bk  = zeros(num,1);
sh  = zeros(num,1);

for n = 1:num
    sttr = trace(strn{n})/dim;
    bk(n) = 0.5*sum(diag(strn{n}).^2);
    sh(n) = 0.5*sum(sum((strn{n}-sttr*eye(dim)).^2));
end