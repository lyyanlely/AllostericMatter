function pr = computePart(vmat,dim)
num = size(vmat,1)/dim;
vmag = zeros(num,size(vmat,2));
for d=1:dim
    vmag = vmag+vmat(d:dim:end,:).^2;
end
pr = 1/num./sum(vmag.^2);
