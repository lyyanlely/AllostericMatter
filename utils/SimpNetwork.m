function [pos, nb] = SimpNetwork(num,dim,z,r)
% this function build networks with different coordination number regions

    ll = round(sqrt(num));
    [pos0,nb0,rs] = ReadNetwork(ll,ll,dim,z,r);
    %L = max(pos0(:,1))-min(pos0(:,1));
    
    shif = repmat([1,0],num,1);
    pos  = pos0;
    
    nr0   = size(nb0,1);
    nr1  = size(rs{1},1);
    nnb  = nr0+nr1;
    nb   = zeros(nnb,2);
    
    for n=1:nnb
        n0 = n;
        if n0<=nr0
            nb1 = nb0(n0,1);
            nb2 = nb0(n0,2);
        else
            n0 = n0-nr0;
            nb1 = rs{1}(n0,1);
            nb2 = rs{1}(n0,2);
        end
        nb(n,:) = [nb1,nb2];
    end