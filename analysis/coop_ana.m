dirc  = '../Allostery/Nov27/';
cname = 'Config';
dtype = '.dat';
dd    = 1;
dt    = 1;
allos = 1;
kweak = 0.0001;
dname = sprintf('%.2f',dd);
tname = sprintf('%.2f',dt);
kname = sprintf('%.4f',kweak);

%npara = 1;
rp    = 20;
beta = 0.01; 

dim = 2;
nx  = 20; %
ny  = 20; %
nlnk = size(bondnb,1); %3*num(j)-2*nx(j);

lx  = length(nx);
npara = length(beta);
pos  = [posx;posy];
%% compute the thermodynamic factors
num  = zeros(lx,1);
nsp  = zeros(lx,1);
zloc = cell(lx,1);  % local coordination number
magn = cell(lx,1);  % magnitude of response
mdisp= cell(lx,1);  % mean displacement
bfac = cell(lx,1);  % b-factor
sbfc = cell(lx,1);  % strain b factor
shre = cell(lx,1);  % shear strain
blke = cell(lx,1);  % bulk strain
strs = cell(lx,1);  % stress
cons = cell(lx,1);  % conservation
flct = cell(lx,1);  % cost of flipping state

for j = 1:lx
    xname  = sprintf('%d',nx(j));
    yname  = sprintf('%d',ny(j));
    num(j) = ny(j)*nx(j);
    nsp(j) = round(num(j)*5/2); 
    nname  = sprintf('%d',nsp(j)); %
    zloc{j} = zeros(npara,num(j));
    magn{j} = zeros(npara,num(j));
    mdisp{j}= zeros(npara,2*num(j));
    bfac{j} = zeros(npara,num(j));
    sbfc{j} = zeros(npara,num(j));
    shre{j} = zeros(npara,num(j));
    blke{j} = zeros(npara,num(j));
    strs{j} = zeros(npara,num(j));
    cons{j} = zeros(npara,nlnk);
    flct{j} = zeros(npara,nlnk);
    nr      = zeros(npara,1);
    for n = 1:npara
        tnam = sprintf('%.4f',beta(n));
        ztemp  = zeros(1,num(j));  % coordination number
        bftemp = zeros(1,num(j));  % b-factor
        sbftmp = zeros(1,num(j));  % strain b-factor
        bktemp = zeros(1,num(j));   % bulk energy
        shtemp = zeros(1,num(j));   % shear energy
        sstemp = zeros(1,num(j));   % stress
        sttemp = zeros(1,nlnk);     % state
        %vortmp = vortmp+vor;  % vorticity
        mdtemp = zeros(1,dim*num(j));
        magtmp = zeros(1,num(j));
        for r = 20 %1:rp
            rname = sprintf('%03d',r);
            filename = [dirc,cname,'_',xname,'_',yname,'_',nname,'_',dname,'_',tname,'_',tnam,'_',kname,'_1_',rname,dtype];  
            if exist(filename,'file')
                s = dir(filename);
                if s.bytes>1
                    cdata = dlmread(filename);
                    lc = size(cdata,1);
                    for l = 91:lc
                        state = cdata(l,:);
                        sttemp = sttemp+state;
                        nb = bondnb(state>0,:);
                        ztemp = ztemp+histc(reshape(nb,1,[]),1:num(j));
                        
                        Mmat  = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
                        [bfc,stb]  = compBfactor(Mmat,pos,bondnb);  % strain b factor
                        if allos == 1  %allosteric cost
                            Qmat   = zeros(dim*num);
                            for i=idd
                                Qmat(i,i) = 1;
                            end
                            Qmat(:,nid) = -Mmat(:,nid);
                            fext = Mmat*pert;
                            disp = linsolve(Qmat,fext);  %kweak>0
                            if any(isnan(disp)|abs(disp)>1e5)
                                disp = pinv(Qmat)*fext;
                            end
                            disp(idd) = pert(idd);
                        else
                            Qmat = zeros(dim*num);
                            for i=[idd(1:2*dim),itt(1:2*dim)]
                                Qmat(i,i) = 1;
                            end
                            Qmat(nidt,nidt) = -Mmat(nidt,nidt);
                            Qmat(idt,nidt)  = -PTmat*Mmat(idt,nidt);
                            Qmat(nidt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)]) = -Mmat(nidt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
                            Qmat(idt,[idd(2*dim+1:4*dim),itt(2*dim+1:4*dim)])  = -PTmat*Mmat(idt,idt)*PTmat([2*dim+1:4*dim,6*dim+1:8*dim],:)';
                            fext = Mmat*(pert+targ);
                            fext(idt) = PTmat*fext(idt);
                            disp = linsolve(Qmat,fext);  %kweak>0
                            if any(isnan(disp)|abs(disp)>1e5)
                                disp = pinv(Qmat)*fext; %
                            end
                            disp([idd(1:2*dim),itt(1:2*dim)]) = 0;
                            disp(idt) = PTmat'*disp(idt);
                            disp = disp+pert+targ;
                        end
                        displ = rmvtsrot(disp,pos); % displacement field with translations and rotations removed
                        stres = diag(diag(state)*Smatrix*displ);
                        sss = cell(num(j),1);
                        for i = 1:num(j)
                            %bid = bondnb(:,1)==i|bondnb(:,2)==i;
                            smm  = Smatrix(:,(i-1)*dim+(1:dim));
                            sss{i} = smm'*stres*smm; %sum(abs(stres(bid)))/2;
                        end
                        [bks,shs] = strainEnergy(sss);
                        [strn,vor] = computeStrain(displ,pos,bondnb);
                        [bk,sh] = strainEnergy(strn);
                        bftemp = bftemp+bfc';  % b-factor
                        sbftmp = sbftmp+stb';  % strain b-factor
                        bktemp = bktemp+bk';   % bulk energy
                        shtemp = shtemp+sh';   % shear energy
                        sstemp = sstemp+bks'+shs';   % stress
                        %vortmp = vortmp+vor;  % vorticity
                        mdtemp = mdtemp+displ';
                        magtmp = magtmp+sqrt(sum(reshape(displ,dim,num).^2));
                        %[cc,~] = histc(sh,xx);
                        %shistt = shistt+cc;
                        nr(n) = nr(n)+1;
                    end
                end
            end
        end
        zloc{j}(n,:) = ztemp/nr(n);
        bfac{j}(n,:) = bftemp/nr(n);
        sbfc{j}(n,:) = sbftmp/nr(n);
        blke{j}(n,:) = bktemp/nr(n);
        shre{j}(n,:) = shtemp/nr(n);
        strs{j}(n,:) = sstemp/nr(n);
        mdisp{j}(n,:)= mdtemp/nr(n);
        magn{j}(n,:) = magtmp/nr(n);
        sttemp = sttemp/nr(n);
        msig   = nsp(j)/nlnk;
        cons{j}(n,:) = sttemp.*log(sttemp/msig)+(1-sttemp).*log((1-sttemp)/(1-msig));
        
    end
    
    % compute the eigenmodes
end

%% compute the soft boundary and define space


