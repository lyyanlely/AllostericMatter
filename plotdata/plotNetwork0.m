function plotNetwork0(pos,nb,xl,flg)
% this function generate plot of network
    if ~exist('flg','var')
    % third parameter does not exist, so default it to something
        flg = zeros(1,size(nb,1))>0;
    end
    %figure; %subplot(1,2,1);
    theta = pi/10;
    x0    = 0.13;
    ccode0= [40/255 150/255 210/255]; %[210/255 100/255 100/255];
    ccode = [210/255 150/255 40/255];
    ccode2= ccode*0.5;
    
    nc = size(nb,1);
    %nr = size(rs,1);
    %pos(:,1) = pos(:,1)-0.5;
    %pos(pos(:,1)<0,1) = pos(pos(:,1)<0,1)+1;
    %plot([.5,.5],[0,1],'--','linewidth',2,'color',[0.85,0.4,0]);
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
            if ~flg(n)
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
            else
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',ccode0,'linewidth',2);
                %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,ccode0,1.5);
            end
            %end
        else
            while dx>.5*xl
                dx = dx-xl;
            end
            while dx<-.5*xl
                dx = dx+xl;
            end
            if dy>0.5
                dy = dy-xl;
            elseif dy<-0.5
                dy = dy+xl;
            end
            %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
            %putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,[120/255 80/255 0],1.5);
            %if (pos(n1,1)-x0)*(pos(n1,1)+dx-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n1,1)+dx-1+x0)<=0 || ...
            %        (pos(n1,2)-x0)*(pos(n1,2)+dy-x0)<0 || (pos(n1,2)-1+x0)*(pos(n1,2)+dy-1+x0)<0
            %    plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
            %    plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'--','color',[0 10/255 240/255],'linewidth',3);
            %else
            if ~flg(n)
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'-','color',ccode,'linewidth',2);
                plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'-','color',ccode,'linewidth',2);
            else
                plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',ccode0,'linewidth',2);
                plot([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],'--','color',ccode0,'linewidth',2);
                %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,ccode0,1.5);
                %putSpring([pos(n2,1) pos(n2,1)-dx],[pos(n2,2) pos(n2,2)-dy],theta,ccode0,1.5);
            end
            %end
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
    %plot(pos(:,1),pos(:,2),'ko','markersize',6,'markerfacecolor',ccode,'color',ccode2,'linewidth',0.5);
    axis equal off;
    %xlim([0,xl]);
    %ylim([0,xl]);
    %set(gca,'xtick',[]);set(gca,'ytick',[]); box on;
    
    %subplot(1,2,2);
    %theta = pi/10;
    %x0    = 0.13;
    %ccode = [210/255 200/255 100/255];
    %ccode2= ccode*0.5;
    
end
