
dirc  = '../Allostery/magn/';
fname = 'singval_fix_';
pname = 'position_fix_';
wname = 'amplitude_fix_';
nname = '1024_';
zname = '4.250_';
dtype = '.dat';

L  = 32;
nb = 2*L;
%ev    = [];
%[~,ix] = findpeaks(sv);
d  = [];
u  = [];
ns = [];
xx = [];
amp= [];
f  = 0;
n = 0;
fracn = sprintf('%.2f_',f);
for x = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] %[0.2,0.4,0.6,0.8] %,1.0,1.4,1.8,2.2,2.6,3.,3.4,3.8]
    xname = sprintf('%.2f_',x);
    for r = 1:100
        flag = 0;
        rname = sprintf('%03d',r);
        %{
        filename = [dirc,wname,nname,xname,fracn,rname,dtype]; %
        if exist(filename,'file')
            s = dir(filename);
                n = n+1;
            if s.bytes>1
                dd = dlmread(filename);
                d  = [d;dd];
                ns = [ns,length(dd)];
                flag = 1;
            end
        end
        %
        filename = [dirc,fname,nname,zname,xname,fracn,rname,dtype];
        if exist(filename,'file')
            dd = dlmread(filename);
            d  = [d;dd];  %(2,:)
            %u  = [u,dd(1,:)];
        end
        %
        filename = [dirc,pname,nname,zname,xname,fracn,rname,dtype];
        if exist(filename,'file')
            data = dlmread(filename);
            %xx  = [xx,data(1,:)]; 
            xx = [xx;data];
        end
        %}
        filename = [dirc,wname,nname,zname,xname,fracn,rname,dtype];
        if exist(filename,'file')
            data = dlmread(filename);
            xx  = [xx,data(1,:)]; 
            amp = [amp,data(2,:)];
        end
        %}
    end
end
%% no lamb dependence
ii=21;
xx1 = reshape(xx,1,[]);
amp1= reshape(amp,1,[]);
amp1= amp1(xx1>0);
xx1 = xx1(xx1>0);
%% compute the distance dependence of magnificaiton
xb = -1/2/nb:1/nb:1+1/2/nb;
[~,ix] = histc(xx1,xb);
aa = zeros(1,nb+1);
va = zeros(1,nb+1);
for i=1:(nb+1)
flg = ix==i;
aa(i)  = mean(log10(amp1(flg&amp1>0)));  %log10
va(i)  = std(log10(amp1(flg&amp1>0)));
end
%sv = [ev;1./ev];
%% compute the size of amplification layer, by half 
xfit = 0:1/2/nb:1;
yfit = spline(0:1/nb:1,aa,xfit);
ap   = (aa(2:end)-aa(1:end-1))*nb/L;
[apo,~] = sort(ap,'descend');
am   = mean(apo(1:2));
am2  = am/2;
yp   = (yfit(2:end)-yfit(1:end-1))*2*nb/L;
id   = find((yp(2:end)-am2).*(am2-yp(1:end-1))>0);
im = id(id>2*nb-8);
i0 = setdiff(id,[im,1]);
if isempty(im) % no open boundary
    im = 2*nb;
end
xm = xfit(im+1);
x0 = mean(xfit(i0+1));
xa = (xm-x0)*L; % size of amplification;
va = std(xfit(i0+1))*L;