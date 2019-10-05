function [pos,nb,fl,fr,nf] = cutNetworkfix1(num,dim,f,r,x0)
% this function build networks with constraints from periodic boundary
    
    lx = round(sqrt(num));
    [pos0,nb0,rs] = ReadNetwork(lx,lx,dim,4.125,r);
    L = max(pos0(:,1))-min(pos0(:,1));
    
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
    idx  = zeros(num,1);
    nc   = num;
    
    nb   = zeros(nr0+nr1,2);
    pos  = zeros(2*num,2);
    pos(1:num,:) = pos0;
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
        dy = (pos(nb2,2)-pos(nb1,2));
        if dx>L/2
            dx = dx-1;
        elseif dx<-L/2
            dx = dx+1;
        end
        if (pos(nb1,1)-x0)*(pos(nb1,1)+dx-x0)>0&&(pos(nb2,1)-x0)*(pos(nb2,1)-dx-x0)>0
            if (dy>L/2||dy<-L/2) && rand(1)<f
                if idx(nb1)>0
                    nb11 = idx(nb1);
                else
                    nc = nc+1;
                    nb11 = nc;
                    idx(nb1)  = nc;
                    pos(nc,:) = pos(nb1,:)-[0,sign(pos(nb1,2)-0.5)];
                end
                if idx(nb2)>0
                    nb22 = idx(nb2);
                else
                    nc = nc+1;
                    nb22 = nc;
                    idx(nb2)  = nc;
                    pos(nc,:) = pos(nb2,:)-[0,sign(pos(nb2,2)-0.5)];
                end
                nnb = nnb+1;
                nb(nnb,:) = [nb1,nb22];
                nnb = nnb+1;
                nb(nnb,:) = [nb2,nb11];
            else
                nnb = nnb+1;
                nb(nnb,:) = [nb1,nb2];
            end
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
    pos = pos(1:nc,:);
    nf   = (num+1):nc;
    fl = setdiff(fl,0);
    fr = setdiff(fr,0);