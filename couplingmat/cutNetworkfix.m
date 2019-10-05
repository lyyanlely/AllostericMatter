function [pos,nb,fl,fr,nf] = cutNetworkfix(num,z,f,r,x0)
% this function build networks with different coordination number regions
    dim = 2;
    lx = round(sqrt(num));
    [pos,nb0,rs] = ReadNetwork(lx,lx,dim,z,r);
    L = max(pos(:,1))-min(pos(:,1));
    
    nfix = round(num*f);
    nf   = zeros(nfix,1);
    
    % divide the system into blocks, so that the randomly fixed particles
    % are homogeneously distributed.
    lb   = max(1,floor(sqrt(nfix)));
    lowb = 0:1/lb:1-0.5/lb;
    uppb = 1/lb:1/lb:1;
    %bset = cell(lb^2,1);
    c0   = floor(nfix/lb^2);
    c1   = mod(nfix,lb^2);
    fexb = randsample(lb^2,c1);
    cc   = 0;
    for ii = 1:lb
        for jj = 1:lb
            bset = find(pos(:,1)>=lowb(ii)&pos(:,1)<uppb(ii)&pos(:,2)>=lowb(jj)&pos(:,2)<uppb(jj));
            if ismember(lb*(ii-1)+jj,fexb)
                nf(cc+1:cc+c0+1) = randsample(bset',c0+1);
                cc = cc+c0+1;
            else
                nf(cc+1:cc+c0) = randsample(bset',c0);
                cc = cc+c0;
            end 
        end
    end
    
    %{
    ndeg = histc(reshape(nb0,1,[]),1:num);
    % remove one connection to the fixed node to keep isostaticity
    for n=1:nfix
        id = sum(ismember(nb0,nf(n)),2);
        neighb = setdiff(reshape(nb0(id>0,:),1,[]),nf(n));
        [~,n0] = sort(ndeg(neighb),'descend');
        for ii=1:min(dim,length(n0))
            if neighb(n0(ii))>nf(n)
                nb0 = setdiff(nb0,[nf(n),neighb(n0(ii))],'rows');
            else
                nb0 = setdiff(nb0,[neighb(n0(ii)),nf(n)],'rows');
            end
        end
    end
    %}
    
    nr0  = size(nb0,1);
    nr1  = size(rs{1},1);
    
    nb   = zeros(nr0+nr1,2);
    fl   = zeros(num,1);
    fr   = zeros(num,1);
    
    nnb  = 0;
    cc   = 0;  % count the number of nodes on the boundary
    for n=1:nr0+nr1
        n0 = n;
        if n0<=nr0
            nb1 = nb0(n0,1);
            nb2 = nb0(n0,2);
        else
            n0 = n0-nr0;
            nb1 = rs{1}(n0,1);
            nb2 = rs{1}(n0,2);
        end
        %ds = sqrt(sum((pos(nb1,:)-pos(nb2,:)).^2));
        dx = (pos(nb2,1)-pos(nb1,1));
        if dx>L/2
            dx = dx-1;
        elseif dx<-L/2
            dx = dx+1;
        end
        if (pos(nb1,1)-x0)*(pos(nb1,1)+dx-x0)>0&&(pos(nb2,1)-x0)*(pos(nb2,1)-dx-x0)>0
            nnb = nnb+1;
            nb(nnb,:) = [nb1,nb2];
        else
            cc = cc+1;
            if dx<0
                fl(cc) = nb1;
                fr(cc) = nb2;
            else
                fl(cc) = nb2;
                fr(cc) = nb1;
            end
        end
    end
    nb = nb(1:nnb,:);
    fl = setdiff(fl,0);
    fr = setdiff(fr,0);