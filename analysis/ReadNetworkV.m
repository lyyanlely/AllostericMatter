function [pos,nb,rs] = ReadNetworkV(num,dim,zz,r)
% find the original network and red springs
%num   = 64;
%r = 1;
%dim   = 2;
z0    = 2*dim;%2*dim-2*dim/num;

dirc  = '../CProgram/network/';
nname = 'network';
rname = sprintf('%03d',r);
dtype = '.dat';
% read the network nb=dn
filename = [dirc,nname,'_',sprintf('%d',num),'_',sprintf('%5.3f',z0),'_',rname,dtype];
d0  = dlmread(filename);
% positions of the particles
pos = mp(d0(:,1:2));
nb  = zeros(dim*num,2);
cc  = 0; % count the connections
for n = 1:num
    for pp = 3:10
        n1 = d0(n,pp)+1;
        if n1<n && n1>0
            cc = cc+1;
            nb(cc,:) = [n1,n];
        end
    end
end

%% add red springs
%z     = 4.0;
nz  = length(zz);
rs  = cell(nz,1);
for i = 1:nz
    z = zz(i);
    filename = [dirc,nname,'_',sprintf('%d',num),'_',sprintf('%5.3f',z),'_',rname,dtype];
    d1    = dlmread(filename);
    rs{i} = zeros(round(num*z/2)-dim*num,2);
    cr    = 0;
    for n = 1:num
        for pp = 3:10
            n1 = d1(n,pp)+1;
            if n1<n && n1>0
                newb = [n1,n];
                if ~ismember(newb,nb,'rows')
                    cr = cr+1;
                    rs{i}(cr,:) = newb;
                end
            end
        end
    end
end