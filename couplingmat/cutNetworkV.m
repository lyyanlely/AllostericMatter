function [pos,nb,fl,fr] = cutNetworkV(num,dim,z,r,x0)
% this function build networks with different coordination number regions

    [pos,nb0,rs] = ReadNetworkV(num,dim,z,r);
    L = max(pos(:,1))-min(pos(:,1));
    
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