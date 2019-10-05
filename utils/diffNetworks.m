function [pos, nb] = diffNetworks(num,dim,z1,z2,x0,r)
% this function build networks with different coordination number regions

    ll = round(sqrt(num));
    [pos0,nb0,rs] = ReadNetwork(ll,ll,dim,[z1,z2],r);
    L = max(pos0(:,1))-min(pos0(:,1));
    
    shif = repmat([1,0],num,1);
    pos0(:,1) = rem(pos0(:,1)+x0,1);
    pos  = [pos0;pos0+shif];
    
    nr0   = size(nb0,1);
    nr1  = size(rs{1},1);
    nr2  = size(rs{2},1);
    nnb  = 2*nr0+nr1+nr2;
    nb   = zeros(nnb,2);
    
    for n=1:nnb
        n0 = n;
        if n0<=nr0+nr1
            if n0<=nr0
                nb1 = nb0(n0,1);
                nb2 = nb0(n0,2);
            else
                n0 = n0-nr0;
                nb1 = rs{1}(n0,1);
                nb2 = rs{1}(n0,2);
            end
            fg1 = 0;  % flag of connection
        else
            n0 = n0-nr0-nr1;
            if n0<=nr0
                nb1 = nb0(n0,1);  % neighbor 
                nb2 = nb0(n0,2);
            else
                n0 = n0-nr0;
                nb1 = rs{2}(n0,1);
                nb2 = rs{2}(n0,2);
            end
            fg1 = 1;  % flag of connection
        end
        %ds = sqrt(sum((pos(nb1,:)-pos(nb2,:)).^2));
        dx = abs(pos(nb1,1)-pos(nb2,1));
        if dx>L/2
            fg2 = 1;
        else
            fg2 = 0;
        end
        nb(n,:) = [nb1+fg1*num,nb2+mod(fg1+fg2,2)*num];
    end