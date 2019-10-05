nf = 30;
dispx0 = dispx;
dispy0 = dispy;
pert0  = pert;
for t = 0:nf
    px = posx+t/nf*dispx0';
    py = posy+t/nf*dispy0';
    scrsz = get(groot,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)])
    plot(px,py,'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',3)
    axis equal off
    hold all
    plot(px(px<1)+xbound,py(px<1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerSize',3)
    
    for i = 1:nsp
        n = occ(i);
        n1 = bondnb(n,1);
        n2 = bondnb(n,2);
        ny1 = ceil(n1/xl);
        nx1 = n1 - (ny1-1)*xl;
        ny2 = ceil(n2/xl);
        nx2 = n2 - (ny2-1)*xl;
        %if ismember(i,find(stresd>0))
        cc = [0.9,0.1,0.2];
        %elseif ismember(i,bondinsp)
        if (nx1==1&&nx2==xl)
            plot([px(n1)+xbound,px(n2)],[py(n1),py(n2)],'linewidth',1,'color',cc);
            plot([px(n1),px(n2)-xbound],[py(n1),py(n2)],'linewidth',1,'color',cc);
        elseif nx2==1&&nx1==xl
            plot([px(n2)+xbound,px(n1)],[py(n2),py(n1)],'linewidth',1,'color',cc);
            plot([px(n2),px(n1)-xbound],[py(n2),py(n1)],'linewidth',1,'color',cc);
        else
            plot([px(n1),px(n2)],[py(n1),py(n2)],'linewidth',1,'color',cc);
        end
    end
    
    quiver(posx(id),posy(id),pert0(idd(1:2:7))',pert0(idd(2:2:8))','linewidth',3,'AutoScale','off','color',[0.6,0.,0.8])
    xlim([-.4,nx-.1])
    ylim([-1,sqrt(3)*ny/2+sqrt(3)])
    f=getframe(gca,rect);
    im(:,:,1,t+1) = rgb2ind(f.cdata,map,'nodither');
    close all
end