%[pos,nb] = diffNetworks(256,2,5.0,4.0,1);
num = 1024;
L = sqrt(num);
dim = 2;
z = 5.0;
r = 1;
upp = 0.5+0.75/L; %0.55;
low = 0.5-0.75/L; %0.45;
%uppb= 0.63;
%lowb= 0.37;
kk = 0.025; %strength of the weak springs
x0 = 0.5;
i0 = 2;
%dd = zeros(20,length(kk));
%oo = zeros(20,length(kk));
%qq = zeros(20,length(kk));
% read the network
[pos,nb,Smat0,bpos]  = trilattice(sqrt(num),sqrt(num),0.);
% make the bonds in [low,upp] softer
flg = bpos(:,1)>low&bpos(:,1)<upp; %&(bpos(:,2)<lowb|bpos(:,2)>uppb);
Smat = Smat0;
Smat(flg,:) = sqrt(kk(1))*Smat0(flg,:);
Mmat = full(Smat'*Smat);
[Vref,Dref] = eig(Mmat);  % eigenvalues and eigenvectors of dynamic matrix
pos(:,2) = mod(pos(:,2)+0.5,1);
[disp,pert,fext] = shearResponse([upp,0.5],[low,0.5],Smat,pos,0); % compute response to a shear dipole
pos(:,2) = mod(pos(:,2)+0.5,1);
d = diag(Dref);
dd = sqrt(d); %frequencies
oo = sum((Vref(1:2:end,:).^2+Vref(2:2:end,:).^2).^2)'; %participation ratio
qq = abs(Vref'*disp)/norm(disp); %overlap between displacement and eigendirections
%{
for i = 2:length(kk)
    Smat(flg,:) = sqrt(kk(i))*Smat0(flg,:);
    Mmat = full(Smat'*Smat);
    [V,D] = eig(Mmat);
    pos(:,2) = mod(pos(:,2)+0.5,1);
    [disp,pert,fext] = shearResponse([upp,0.5],[low,0.5],Smat,pos,2);
    pos(:,2) = mod(pos(:,2)+0.5,1);
    d = diag(D);
    dd(:,i) = sqrt(d(i0+(1:20)));
    qq(:,i) = abs(V(:,i0+(1:20))'*disp)/norm(disp);
    [~,id]  = max(abs(V(:,i0+(1:20))'*Vref),[],2);
    for j = 1:20;
        if id(j)-i0<21
            oo(id(j)-i0,i) = sqrt(d(i0+j));
        end
    end
    oo(oo(:,i)==0,i)=nan;
end
%}
%plotNetwork0(pos,nb,1,flg);
%plot(bpos(:,1),bpos(:,2),'bo');
% compute point response
%%
tic;
[~,es,~] = shearResponse([upp,0.5],[low,0.5],Smat,pos,0);  % shear energy to dipole at allosteric site
[~,et,~] = shearResponse([upp,0.5],[low,0.5],Smat,pos,1);  % shear energy to dipole at active site
[~,est,~] = shearResponse([upp,0.5],[low,0.5],Smat,pos,2); % shear energy to dipole at both sites
pcost = es+et-est;
pcost0 = pcost;   % cooperative energy of a reference configuration
toc;
nnb = size(nb,1);
flptmp = zeros(nnb,1);
tic;
for ii=1:nnb
    Smat1 = Smat;
    if flg(ii)
        Smat1(ii,:) = Smat1(ii,:)/sqrt(kk);
    else
        Smat1(ii,:) = sqrt(kk)*Smat1(ii,:);
    end
    %pMmat = Smat1'*Smat1;
    [~,es,~] = shearResponse([upp,0.5],[low,0.5],Smat1,pos,0);
    [~,et,~] = shearResponse([upp,0.5],[low,0.5],Smat1,pos,1);
    [~,est,~] = shearResponse([upp,0.5],[low,0.5],Smat1,pos,2);
    pcost = es+et-est;
    flptmp(ii) = flptmp(ii)+pcost-pcost0;  % energy difference
end
toc;
%%
fcst  = zeros(num,1);  % flipping cost on each node
for i=1:num   % average link values to each node
    ffg = ismember(nb(:,1),i)|ismember(nb(:,2),i);
    fcst(i) = mean(flip(ffg));
end
%%
[bf,stb]  = compBfactor(Mmat,pos',nb);

Smatrix = Smat0;
Smatrix(flg,:) = kk*Smatrix(flg,:);
stres = diag(Smatrix*disp);
sss = cell(num,1);
for i = 1:num
    smm  = Smat(:,(i-1)*dim+(1:dim));
    sss{i} = smm'*stres*smm; %sum(abs(stres(bid)))/2;
end
[bks,shs] = strainEnergy(sss);
[strn,vor] = computeStrain(disp,pos',nb);
[bk,sh] = strainEnergy(strn);
magtmp = sqrt(sum(reshape(disp,dim,num).^2))';
%%
id = union(ceil(find(pert~=0)/2),[]);
quiver(pos(:,1),pos(:,2),disp(1:2:end),disp(2:2:end),'linewidth',.5,'AutoScaleFactor',3,'color',[0.1,0.1,0.1])
hold all;
quiver(pos(id,1),pos(id,2),disp(2*id-1),disp(2*id),'linewidth',2,'AutoScaleFactor',.3,'color',[0.6,0.,0.8])
%plot([0.,1.],[upp,upp],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([0.,1.],[low,low],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([1.,1.],[0,1],'--','linewidth',2,'color',[0.85,0.4,0]);