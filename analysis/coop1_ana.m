
minp = xl/num;
neig = ceil(num*dim/20);
nconf = size(config,1);
spf  = nsp/nnb;
nend = 34;  % number of configurations each sequence
ns  = nconf/nend;  % number of sequences
eps = 1e-15;
%
sid = [];
for i=0:(xl-np)
    sid = [sid;(1:np)+i];
end
sid = [sid;it];
nsid = size(sid,1);

pos = [posx;posy];
%}

s = rng('shuffle');
msig  = zeros(ns,nnb);  % mean occupancy
flip  = zeros(ns,nnb);  % cost of flipping
mdisp = zeros(ns,num*dim);
magn  = zeros(ns,num);
cons  = zeros(ns,num);  % conservation on each node
fcst  = zeros(ns,num);  % flipping cost on each node
eshr  = zeros(ns,num);  % shear pseudoenergy
sbfc  = zeros(ns,num);  % strain B-factor
sencost = zeros(ns,nsid);
senmagn = zeros(ns,nsid*num);
dd = zeros(ns,neig*nend);
pr = zeros(ns,neig*nend);
qq = zeros(ns,neig*nend);
rd = zeros(ns,neig*nend);  % random configurations
rp = zeros(ns,neig*nend);
rq = zeros(ns,neig*nend);
for r=1:ns
    tic;
    %% conservation
    msig(r,:) = mean(config((r-1)*nend+1:r*nend,:));
    msigc   = min(max(msig(r,:),eps),1-eps);
    conserv = msigc.*log(msigc/spf)+(1-msigc).*log((1-msigc)/(1-spf));  % conservation on each link
    for i=1:nend
        state = config((r-1)*nend+i,:);
        %% response
        Mmat = kweak*MmatW+Smatrix'*diag(state)*Smatrix;
        disp = computeResp(Mmat,pos,allos,pert,idd,Pmat);
        nd   = norm(disp);
        disp = disp'/nd;
        mdisp(r,:) = mdisp(r,:)+disp/nend;
        mm = zeros(1,num);
        for d = 1:dim
            mm = mm+disp(d:dim:end).^2;
        end
        magn(r,:) = magn(r,:)+sqrt(mm)/nend*sqrt(num);
        %% pseudo-energy
        [strn,vor] = computeStrain(disp',pos,bondnb);
        [bk,sh] = strainEnergy(strn);
        eshr(r,:) = eshr(r,:)+sh'/nend;
        %% omega, participation, overlap
        [V,d]=eigs(Mmat+eps*eye(num*dim),neig,'sm');
        pr(r,(i-1)*neig+1:i*neig) = computePart(V,2);
        dd(r,(i-1)*neig+1:i*neig) = diag(d)';
        qq(r,(i-1)*neig+1:i*neig) = (disp*V).^2;
        %% B-factor
        if i==nend
            %tic;
            [bf,stb]  = compBfactor(Mmat,pos,bondnb);
            %toc;
            sbfc(r,:) = stb';
        end
        %% sensitivity test
        %
        if i==nend
            %tic;
            for ii=1:nsid
                [pp,tt,ids,mats] = prepareTest(pos,sid(ii,:),it);
                [pcost,displ] = comppcost(Mmat,pp,tt,allos,pos,ids,mats);
                displ = displ'/nd;
                sencost(r,ii) = pcost;
                mm = zeros(1,num);
                for d = 1:dim
                    mm = mm+displ(d:dim:end).^2;
                end
                senmagn(r,(ii-1)*num+1:ii*num) = sqrt(mm)*sqrt(num);
            end
            %toc;
        end
        %}
        %% flipping cost
        if i==nend
            %tic;
            [pcost,~] = comppcost(Mmat,pert,targ,allos,pos,{idd,nid,itt,nit,idt,nidt},{Pmat,Tmat,PTmat});
            pcost0 = pcost;
            for n=1:nnb
                state1 = state;
                state1(n) = mod(state(n)+1,2);
                Mmat = kweak*MmatW+Smatrix'*diag(state1)*Smatrix;
                [pcost,~] = comppcost(Mmat,pert,targ,allos,pos,{idd,nid,itt,nit,idt,nidt},{Pmat,Tmat,PTmat});
                flip(r,n) = pcost0-pcost;
            end
            %toc;
        end
        %% random configuration
        occ = zeros(nnb,1);
        idx = randsample(nnb,nsp);
        occ(idx) = 1;
        Mmat = kweak*MmatW+Smatrix'*diag(occ)*Smatrix;
        [V,d]=eigs(Mmat+eps*eye(num*dim),neig,'sm');
        rp(r,(i-1)*neig+1:i*neig) = computePart(V,2);
        rd(r,(i-1)*neig+1:i*neig) = diag(d)';
        rq(r,(i-1)*neig+1:i*neig) = (disp*V).^2;
    end
    % from bond to site
    for n=1:num   % average link values to each node
        flg = ismember(bondnb(:,1),n)|ismember(bondnb(:,2),n);
        cons(r,n) = mean(conserv(flg));
        fcst(r,n) = mean(flip(r,flg));
    end
    toc;
end
