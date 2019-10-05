function [posl,nb,fl,fr] = cutSlab(num,dim,z,xn,r)
% this function cut the square network into slabs
    %% read network
    [pos,nb0,rs] = ReadNetwork(num,dim,z,r);
    L = 1.; % max(pos(:,1))-min(pos(:,1));
 %% size
    nl = max(1,round(num^(1/dim)/xn));
    dx = L/nl;
    x0 = 0:dx:(L-dx);  % position of boundaries.
    x1 = dx:dx:L;      % upper boundaries
    nr0  = size(nb0,1);
    nr1  = size(rs{1},1);
    
    %% return variables
    posl = cell(nl,1);
    nb   = cell(nl,1);
    fl  = cell(nl,1);
    fr  = cell(nl,1);
    
%% cut slabs
    for i=1:nl
        posl{i} = pos(pos(:,1)>x0(i)&pos(:,1)<=x1(i),:);
        nb{i}   = zeros(nr0+nr1,2);
        fl{i}   = zeros(size(posl{i},1),1);
        fr{i}   = zeros(size(posl{i},1),1);
    end
    
    nnb  = zeros(nl,1);
    cc   = zeros(nl,1);  % count the number of nodes on the boundary
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
            dx = dx-L;
        elseif dx<-L/2
            dx = dx+L;
        end
        %% find the slab
        i1 = find(pos(nb1,1)>x0&pos(nb1,1)<=x1);
        i2 = find(pos(nb2,1)>x0&pos(nb2,1)<=x1);
        fg = (pos(nb1,1)-x0).*(pos(nb1,1)+dx-x0)<0|(pos(nb2,1)-x0).*(pos(nb2,1)-dx-x0)<0; 
        n1 = find(ismember(posl{i1},pos(nb1,:),'rows'));
        n2 = find(ismember(posl{i2},pos(nb2,:),'rows'));
        if all(fg==0) %i1==i2(one slab) % not a bond through the boundary
            nnb(i1) = nnb(i1)+1;
            nb{i1}(nnb(i1),:) = [n1,n2];
        else
            cc(i1) = cc(i1)+1;
            cc(i2) = cc(i2)+1;
            if dx<0
                fl{i1}(cc(i1)) = n1;
                fr{i2}(cc(i2)) = n2;
            else
                fl{i2}(cc(i2)) = n2;
                fr{i1}(cc(i1)) = n1;
            end
        end
    end
    for i=1:nl
        nb{i} = nb{i}(1:nnb(i),:);
        fl{i} = setdiff(fl{i},0);
        fr{i} = setdiff(fr{i},0);
    end