dirc  = '../Allostery/3D/';
dpnam = 'MDisp3D';
bevon = 'Evobond3D';
sfeat = 'SiteFeat3D';
specn = 'Spectrum3D';
kname = sprintf('%.4f',kweak);
aname = sprintf('%d',allos);
xname  = sprintf('%d',xl);
yname  = sprintf('%d',yl);
zname  = sprintf('%d',zl);
nname  = sprintf('%d',nsp);
tnam = sprintf('%.4f',beta);

nr = 10;

msig  = zeros(nr,nnb);  % mean occupancy
flip  = zeros(nr,nnb);  % cost of flipping
mdisp = zeros(nr,num*dim);
magn  = zeros(nr,num);
cons  = zeros(nr,num);  % conservation on each node
fcst  = zeros(nr,num);  % flipping cost on each node
eshr  = zeros(nr,num);  % shear pseudoenergy
sbfc  = zeros(nr,num);  % strain B-factor
dd = zeros(nr,neig*nend);
pr = zeros(nr,neig*nend);
qq = zeros(nr,neig*nend);
rd = zeros(nr,neig*nend);  % random configurations
rp = zeros(nr,neig*nend);
rq = zeros(nr,neig*nend);

xid = [];
for r = 1:nr
    rname = sprintf('%03d',r);
    mdpname = [dirc,dpnam '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
    bevname = [dirc,bevon '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
    sftname = [dirc,sfeat '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
    spcname = [dirc,specn '_' aname '_' xname '_' yname '_' zname '_' nname '_' tnam '_' kname '_' rname '.dat'];
    if exist(mdpname,'file')
        xid = [xid,r];
        mdisp(r,:) = dlmread(mdpname);
        bev = dlmread(bevname);
        msig(r,:)  = bev(1,:);
        flip(r,:)  = bev(2,:);
        sft = dlmread(sftname);
        magn(r,:)  = sft(1,:);
        cons(r,:)  = sft(2,:);
        fcst(r,:)  = sft(3,:);
        eshr(r,:)  = sft(4,:);
        sbfc(r,:)  = sft(5,:);
        spc = dlmread(spcname);
        dd(r,:) = spc(1,:);
        pr(r,:) = spc(2,:);
        qq(r,:) = spc(3,:);
        rd(r,:) = spc(4,:);
        rp(r,:) = spc(5,:);
        rq(r,:) = spc(6,:);
    end
end
