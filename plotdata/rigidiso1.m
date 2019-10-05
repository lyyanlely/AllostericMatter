function [pos,nb,fl,fr] = rigidiso1(num,dim,z1,z2,x0,r)
% this function build networks with different coordination number regions

    %z2 = 4;
    ll = round(sqrt(num));
    [pos0,nb0,rs] = ReadNetwork(ll,ll,dim,[z1,z2],r);
    L = max(pos0(:,1))-min(pos0(:,1));
    
    shfx = repmat([1,0],num,1);
    shfy = repmat([0,1],num,1);
    shif = repmat([1,1],num,1);
    pos0(:,1) = rem(pos0(:,1)+x0,1);
    pos  = [pos0;pos0+shfx;pos0+shfy;pos0+shif];
    
    nr0   = size(nb0,1);
    nr1  = size(rs{1},1);
    nr2  = size(rs{2},1);
    nnb  = 4*nr0+3*nr1+nr2;
    nb   = zeros(nnb,2);
    fl   = zeros(num,1);
    fr   = zeros(num,1);

    nbn  = 0;
    cc   = 0;
    for n=1:nnb
        if n<=nr0+nr1 || n>2*nr0+nr1+nr2
            if n<=nr0+nr1
                n0 = n;
                fg1 = 0;  % flag of connection
            else
                n0 = mod(n-nr0-nr2-1,nr0+nr1)+1;
                if n<=3*nr0+2*nr1+nr2
                    fg1 = 2;  % flag of connection
                else
                    fg1 = 3;
                end
            end
            if n0<=nr0
                nb1 = nb0(n0,1);
                nb2 = nb0(n0,2);
            else
                n0 = n0-nr0;
                nb1 = rs{1}(n0,1);
                nb2 = rs{1}(n0,2);
            end
        else
            n0 = n-nr0-nr1;
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
        dx = (pos(nb1,1)-pos(nb2,1));
        if dx>L/2
            fg2 = 1;
        elseif dx<-L/2
            fg2 = -1;
        else
            fg2 = 0;
        end
        dy = abs(pos(nb1,2)-pos(nb2,2));
        if dy>L/2
            fg3 = fg2+2; % in the second row
        else
            fg3 = fg2;  % in the first row
        end
        %
        %if fg2 && ((~fg1 && pos(nb2,1)>0.5)||(fg1 && pos(nb2,1)<0.5))
        if fg2<0 && mod(fg1,2)==0
            cc = cc+1;
            fl(cc) = nb1+fg1*num;
            fr(cc) = nb2+mod(fg1+fg3+4,4)*num;
        elseif fg2>0 && mod(fg1,2)==1
            cc = cc+1;
            fl(cc) = nb2+mod(fg1+fg3+4,4)*num;
            fr(cc) = nb1+fg1*num;
        else
        %else
        %}
            nbn=nbn+1;
            nb(nbn,:) = [nb1+fg1*num,nb2+mod(fg1+fg3+4,4)*num];
        end
    end
    nb = nb(1:nbn,:);
    fl = setdiff(fl,0);
    fr = setdiff(fr,0);