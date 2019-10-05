function make_two_networks_fig

    scrsz = get(0,'ScreenSize');
    fig1 = figure('Position',[scrsz(3)/4 scrsz(4)/5 scrsz(4)/1.5 scrsz(4)/5.0]);

    data = load('network_256_4.000_000.dat');
    %data = data(2:end,:);
    
    N = length(data);
    theta = pi/10;
    
    close all;
    
    figure; subplot(1,2,1);
    hold on;
    for i=1:N-1
        j=3;
        neb = data(i,j)+1;
        while ( neb ~= 0 )
            if ( neb > i ) %no repeating please
                x1 = data(i,1);
                x2 = data(neb,1);
                y1 = data(i,2);
                y2 = data(neb,2);
                dx = x2 - x1;
                dy = y2 - y1;
                flag = 1;
                if (dy > 0.5)
                    flag = 0;
                    dy = dy - 1.0;
                elseif (dy < -0.5)
                    flag = 0;
                    dy = dy + 1.0;
                end
                if (dx > 0.5)
                    flag = 0;
                    dx = dx - 1.0;
                elseif (dx < -0.5)
                    flag = 0;
                    dx = dx + 1.0;
                end
                putSpring([x1 x1+dx],[y1 y1+dy],theta,[120/255 80/255 0],2.5);
                if (flag==0)
                    putSpring([x2 x2-dx],[y2 y2-dy],theta,[120/255 80/255 0],2.5);
                end
            end
            j = j+1;
            neb = data(i,j)+1;
        end
    end
    %[255/256 153/256 0],'color',[255/256 153/256 0]);
    plot(data(:,1),data(:,2),'ko','markersize',10,'markerfacecolor',[255/255 170/255 0],'color',[120/255 80/255 0],'linewidth',0.5);
    axis equal;
    xlim([0.2 0.7]);
    ylim([0.2 0.7]);
    set(gca,'xtick',[]);set(gca,'ytick',[]); box on;


    %{
    data = load('hex.dat');
    data = data(2:end,:);
    
    N = length(data);
    theta = pi/10;
    %}
    x = data(:,1);
    y = data(:,2);
    subplot(1,2,2);
    hold on;
    for i=1:N-1
                x1 = data(i,1);
                y1 = data(i,2);
        if x1>=0.2 && x1<=0.7 && y1>=0.2 && y1<=0.7
        j=3;
        nb  = data(i,3:end)+1;
        neb = data(i,j)+1;
        while ( neb ~= 0 )
            if ( neb > i ) %no repeating please
                x2 = data(neb,1);
                y2 = data(neb,2);
                dx = x2 - x1;
                dy = y2 - y1;
                flag = 1;
                if (dy > 0.5)
                    flag = 0;
                    dy = dy - 1.0;
                elseif (dy < -0.5)
                    flag = 0;
                    dy = dy + 1.0;
                end
                if (dx > 0.5)
                    flag = 0;
                    dx = dx - 1.0;
                elseif (dx < -0.5)
                    flag = 0;
                    dx = dx + 1.0;
                end
                putSpring([x1 x1+dx],[y1 y1+dy],theta,[240/255 220/255 190/255],2.5);
                if (flag==0)
                    putSpring([x2 x2-dx],[y2 y2-dy],theta,[240/255 220/255 190/255],2.5);
                end
            end
            j = j+1;
            neb = data(i,j)+1;
        end
        dist = sqrt((x-x1).^2+(y-y1).^2);
        [~,ix] = sort(dist);
        cc = 1;
        ii = 2;
        while cc<=6
            if ~ismember(ix(ii),nb)
                x2 = data(ix(ii),1);
                y2 = data(ix(ii),2);
                dx = x2 - x1;
                dy = y2 - y1;
                flag = 1;
                if (dy > 0.5)
                    flag = 0;
                    dy = dy - 1.0;
                elseif (dy < -0.5)
                    flag = 0;
                    dy = dy + 1.0;
                end
                if (dx > 0.5)
                    flag = 0;
                    dx = dx - 1.0;
                elseif (dx < -0.5)
                    flag = 0;
                    dx = dx + 1.0;
                end
                plot([x1 x1+dx],[y1 y1+dy],'color',[190/255 220/255 240/255],'linewidth',1.5);
                if (flag==0)
                    plot([x2 x2-dx],[y2 y2-dy],'color',[190/255 220/255 240/255],'linewidth',1.5);
                end
            cc = cc+1;
            end
            ii = ii+1;
        end
        end
    end
    %[255/256 153/256 0],'color',[255/256 153/256 0]);
    plot(data(:,1),data(:,2),'ko','markersize',25,'markerfacecolor',[255/255 170/255 0],'color',[120/255 80/255 0],'linewidth',0.5);
    
    axis equal;
    xlim([0.35 0.6]);
    ylim([0.41 0.66]);

    %}
%       axis([-0.0025 1.0025 -0.0025 1.0025]); box on;
%     title (['file: ',filename, ', ',num2str(n),' particles']);
     set(gca,'xtick',[]);set(gca,'ytick',[]); box on;
     hold off;
     
     h1 = subplot(1,2,1);
     h2 = subplot(1,2,2);
     
     set(h1,'Position',[0.05 0.05 0.425 0.9]);
     set(h2,'Position',[0.55 0.05 0.425 0.9]);
     
end

function plotParticle(r,x,y,color,fill)
    rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1,1],'EdgeColor',color,'LineWidth',1.5,'FaceColor',fill);
end
   


function putSpring(rx,ry,theta,color,linewidth)
    dx = rx(2) - rx(1);
    dy = ry(2) - ry(1);
    r = sqrt(dx*dx + dy*dy);
    dx = dx/r; dy = dy/r;
    foo = r/6;
    
    h = foo*tan(theta);
    sin_alpha = h/sqrt(h^2 + 0.25*foo^2);
    cos_alpha = 0.5*foo/sqrt(h^2 + 0.25*foo^2);
    a = sqrt(foo^2 + h^2);
    b = sqrt(foo^2 + 4*h^2);
    
    x(1) = rx(1); y(1) = ry(1);
    x(7) = rx(2); y(7) = ry(2);
    x(2) = rx(1) + a*(dx*cos(theta)-dy*sin(theta));
    y(2) = ry(1) + a*(dy*cos(theta)+dx*sin(theta));
    x(6) = rx(2) + a*(dx*cos(pi-theta)-dy*sin(pi-theta));
    y(6) = ry(2) + a*(dy*cos(pi-theta)+dx*sin(pi-theta));
  
    
    x(3) = x(2) + b*(dx*cos_alpha+dy*sin_alpha);
    y(3) = y(2) + b*(dy*cos_alpha-dx*sin_alpha);
    
    x(4) = x(3) + b*(dx*cos_alpha-dy*sin_alpha);
    y(4) = y(3) + b*(dy*cos_alpha+dx*sin_alpha);
    
    x(5) = x(4) + b*(dx*cos_alpha+dy*sin_alpha);
    y(5) = y(4) + b*(dy*cos_alpha-dx*sin_alpha);
    
%     plot([x(1) x(2)],[y(1) y(2)],'color',color,'linewidth',linewidth);
%     plot([x(2) x(3)],[y(2) y(3)],'color',color,'linewidth',linewidth);
%     plot([x(3) x(4)],[y(3) y(4)],'color',color,'linewidth',linewidth);
%     plot([x(4) x(5)],[y(4) y(5)],'color',color,'linewidth',linewidth);
%     plot([x(5) x(6)],[y(5) y(6)],'color',color,'linewidth',linewidth);
%     plot([x(6) x(7)],[y(6) y(7)],'color',color,'linewidth',linewidth);

    plot(x,y,'--','color',color,'linewidth',linewidth);
    
end
