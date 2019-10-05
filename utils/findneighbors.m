function nb=findneighbors(pos,th,pflg)
num = size(pos,1);
dim = size(pos,2);
xbound = ceil(max(pos,1));
nb = zeros(12*num,2);
nnb = 0;
for i = 1:num-1
    for j = i+1:num
        dx = zeros(1,dim);
        for d = 1:dim
            dx(d) = pos(i,d)-pos(j,d);
            if pflg
                if dx(d)>xbound/2
                    dx(d) = dx(d)-xbound;
                elseif dx(d)<-xbound/2
                    dx(d) = dx(d)+xbound;
                end
            end
        end
        dr = norm(dx);
        if dr<th
            nnb = nnb+1;
            nb(nnb,1) = i;
            nb(nnb,2) = j;
        end
    end
end
nb = nb(1:nnb,:);        