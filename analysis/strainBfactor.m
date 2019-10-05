function stb = strainBfactor(Mmat,pos,nb,opt)
% this function compute the strain b-factor
eps = 1e-8;  % minimal eigen value contributed
dim = size(pos,1);
num = size(pos,2);

stb = zeros(num,1);

[V,D] = eig(Mmat);

lamb = diag(D);

tic;
for i = 1:dim*num
    if lamb(i)>eps
        disp = V(:,i);
        if nargin>3
            [strn,~] = computeStrain(disp,pos,nb,opt);
        else
            [strn,~] = computeStrain(disp,pos,nb);
        end
        [bk,sh] = strainEnergy(strn);
        stb = stb+2*(bk+sh)/lamb(i);
    end
end
toc;