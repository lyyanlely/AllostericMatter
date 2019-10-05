%config =dlmread('../Allostery/Apr6/Config_12_12_270_1.00_1.00_0.0001_9_001.dat');
%pert(idd([1,2,7,8])) = ptmp(idd([1,2,7,8]));
i = 72;
stat = pconf(i,:); 
Mmat = kweak*MmatW+Smatrix'*diag(stat)*Smatrix;
%stat = stat(1:nnb);
occ = find(stat>0);
vac = find(stat==0);
%posx = rem(posx+4.,nx);
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)])
plot(posx,posy,'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',3)
axis equal off
hold all
plot(posx(posx<1)+xbound,posy(posx<1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerSize',3)
%[stresd,bondinsp,rigidset,pivot] = allo_rigidity(yl,xl,nsp,[posx;posy]',bondnb(occ,:));
%
for i = 1:nsp
    n = occ(i);
    n1 = bondnb(n,1);
    n2 = bondnb(n,2);
    ny1 = ceil(n1/xl);
    nx1 = n1 - (ny1-1)*xl;
    ny2 = ceil(n2/xl);
    nx2 = n2 - (ny2-1)*xl;
    if ~flg(n) %ismember(i,find(stresd>0))
        cc = [0.9,0.1,0.2];
    %elseif ismember(i,bondinsp)
    %    cc = [0.2,0.9,0.1];
    else
        cc = [0.1,0.2,0.9];
    end
    if (nx1==1&&nx2==xl)
        plot([posx(n1)+xbound,posx(n2)],[posy(n1),posy(n2)],'linewidth',1,'color',cc);
        plot([posx(n1),posx(n2)-xbound],[posy(n1),posy(n2)],'linewidth',1,'color',cc);
    elseif nx2==1&&nx1==xl
        plot([posx(n2)+xbound,posx(n1)],[posy(n2),posy(n1)],'linewidth',1,'color',cc);
        plot([posx(n2),posx(n1)-xbound],[posy(n2),posy(n1)],'linewidth',1,'color',cc);
    else
        plot([posx(n1),posx(n2)],[posy(n1),posy(n2)],'linewidth',1,'color',cc);
    end
end

quiver(posx(id),posy(id),pert(idd(1:2:7))',pert(idd(2:2:8))','linewidth',3,'AutoScale','off','color',[0.6,0.,0.8])
%disp   = disp+pconst(1)*pfx+pconst(2)*pfy+pconst(3)*pfr;
[pest,displ] = comppcost(Mmat,pert,targ,allos,[posx;posy],{idd,nid,itt,nit,idt,nidt},{Pmat,Tmat,PTmat});
dispx = displ(1:2:end);
dispy = displ(2:2:end);
%}
quiver(posx(tid)',posy(tid)',dispx(tid),dispy(tid),'linewidth',1.5,'AutoScale','off','color',[0.1,0.1,0.1])

xlim([-.4,nx-.1])
ylim([-1,sqrt(3)*ny/2+sqrt(3)])