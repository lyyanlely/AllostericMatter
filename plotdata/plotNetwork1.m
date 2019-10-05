function plotNetwork1(pos,nb,nr)
% this function generate plot of network
    figure; subplot(1,2,1);
    theta = pi/10;
    x0    = 0.13;
    ccode = [210/255 200/255 100/255];
    ccode2= ccode*0.5;
    
    nc = size(nb,1);
    %nr = size(rs,1);
    
    hold on;
    for n=1:nc
        n1 = nb(n,1);
        n2 = nb(n,2);
        dx = pos(n2,1)-pos(n1,1);
        dy = pos(n2,2)-pos(n1,2);
        if abs(dx)<0.5 && abs(dy)<0.5
            %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
            if (pos(n1,1)-x0)*(pos(n2,1)-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n2,1)-1+x0)<=0 || ...
                    (pos(n1,2)-x0)*(pos(n2,2)-x0)<0 || (pos(n1,2)-1+x0)*(pos(n2,2)-1+x0)<0
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
            else
            plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
            end
        else
            if dx>0.5
                dx = dx-1;
            elseif dx<-0.5
                dx = dx+1;
            end
            if dy>0.5
                dy = dy-1;
            elseif dy<-0.5
                dy = dy+1;
            end
            %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
            %putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,[120/255 80/255 0],1.5);
            if (pos(n1,1)-x0)*(pos(n1,1)+dx-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n1,1)+dx-1+x0)<=0 || ...
                    (pos(n1,2)-x0)*(pos(n1,2)+dy-x0)<0 || (pos(n1,2)-1+x0)*(pos(n1,2)+dy-1+x0)<0
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
                plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'--','color',[0 10/255 240/255],'linewidth',3);
            else
            plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
            plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'-','color',ccode,'linewidth',2);
            end
        end
    end
    % red springs
    %{
    for n=1:nr
        n1 = rs(n,1);
        n2 = rs(n,2);
        dx = pos(n2,1)-pos(n1,1);
        dy = pos(n2,2)-pos(n1,2);
        if abs(dx)<0.5 && abs(dy)<0.5
            putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[240/255 10/255 0],3);
        else
            if dx>0.5
                dx = dx-1;
            elseif dx<-0.5
                dx = dx+1;
            end
            if dy>0.5
                dy = dy-1;
            elseif dy<-0.5
                dy = dy+1;
            end
            putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[240/255 10/255 0],3);
            putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,[240/255 10/255 0],3);
        end
    end
    %}
    plot(pos(:,1),pos(:,2),'ko','markersize',6,'markerfacecolor',ccode,'color',ccode2,'linewidth',0.5);
    axis equal;
    xlim([0,1.]);
    ylim([0,1.]);
    set(gca,'xtick',[]);set(gca,'ytick',[]); box on;
    
    subplot(1,2,2);
    %theta = pi/10;
    %x0    = 0.13;
    %ccode = [210/255 200/255 100/255];
    %ccode2= ccode*0.5;
    
    nc = size(nb,1);
    %nr = size(rs,1);
    
    hold on;
    for n=1:nc
        n1 = nb(n,1);
        n2 = nb(n,2);
        dx = pos(n2,1)-pos(n1,1);
        dy = pos(n2,2)-pos(n1,2);
        if abs(dx)<0.5 && abs(dy)<0.5
            %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
            %if (pos(n1,1)-x0)*(pos(n2,1)-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n2,1)-1+x0)<=0 || ...
            %        (pos(n1,2)-x0)*(pos(n2,2)-x0)<0 || (pos(n1,2)-1+x0)*(pos(n2,2)-1+x0)<0
            %    plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
            %else
            plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
            %end
        else
            if dx>0.5
                dx = dx-1;
            elseif dx<-0.5
                dx = dx+1;
            end
            if dy>0.5
                dy = dy-1;
            elseif dy<-0.5
                dy = dy+1;
            end
            %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
            %putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,[120/255 80/255 0],1.5);
            %if (pos(n1,1)-x0)*(pos(n1,1)+dx-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n1,1)+dx-1+x0)<=0 || ...
            %        (pos(n1,2)-x0)*(pos(n1,2)+dy-x0)<0 || (pos(n1,2)-1+x0)*(pos(n1,2)+dy-1+x0)<0
            %    plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
            %    plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'--','color',[0 10/255 240/255],'linewidth',3);
            %else
            plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
            plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'-','color',ccode,'linewidth',2);
            %end
        end
    end
    % red springs
    %
    nl = round(sqrt(nr));
    nl2= nl^2;
    cc = zeros(nl2,2);
    for n=1:nl2
        ny = floor((n-1)/nl);
        nx = n-nl*ny-1;
        cc(n,1) = (nx+0.5)/nl;
        cc(n,2) = (ny+0.4)/nl;
    end
    for n=1:nl2
        n1   = nearestNeighbor(pos,cc(n,:));
        n2   = [];
        pset = setdiff(pos,pos(n1,:),'rows');
        while isempty(n2)
            nn = nearestNeighbor(pset,cc(n,:));
            ns = find(ismember(pos,pset(nn,:),'rows'));
            if ~ismember([n1,ns],nb,'rows') && ~ismember([ns,n1],nb,'rows')
                n2 = ns;
            end
            pset = setdiff(pset,pset(nn,:),'rows');
        end
        dx = pos(n2,1)-pos(n1,1);
        dy = pos(n2,2)-pos(n1,2);
        if abs(dx)<0.5 && abs(dy)<0.5
            putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[240/255 10/255 0],3);
        else
            if dx>0.5
                dx = dx-1;
            elseif dx<-0.5
                dx = dx+1;
            end
            if dy>0.5
                dy = dy-1;
            elseif dy<-0.5
                dy = dy+1;
            end
            putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[240/255 10/255 0],3);
            putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,[240/255 10/255 0],3);
        end
    end
    %}
    plot(pos(:,1),pos(:,2),'ko','markersize',6,'markerfacecolor',ccode,'color',ccode2,'linewidth',0.5);
    axis equal;
    xlim([0,1.]);
    ylim([0,1.]);
    set(gca,'xtick',[]);set(gca,'ytick',[]); box on;
end
