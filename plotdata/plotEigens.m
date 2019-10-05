ff=3;
ss=1;

%stat= config(end,:,2);
occ = find(stat>0);
vac = find(stat==0);
plot(posx,posy,'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',16)
axis equal off
hold all
plot(posx(posx<1)+xbound,posy(posx<1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerSize',16)
[stresd,bondinsp,rigidset,pivot] = allo_rigidity(yl,xl,nsp,[posx;posy]',bondnb(occ,:));
Nmat = Smatrix(stat'>0,:)*Smatrix(stat'>0,:)';
Mmat = Smatrix'*diag(stat)*Smatrix;
[VV,~] = eig(full(Nmat));
Vu = VV(:,ss);
for i = 1:nsp
    n = occ(i);
    n1 = bondnb(n,1);
    n2 = bondnb(n,2);
    ny1 = ceil(n1/xl);
    nx1 = n1 - (ny1-1)*xl;
    ny2 = ceil(n2/xl);
    nx2 = n2 - (ny2-1)*xl;
    if Vu(i)>0
        cc = [0.7,0.6,0.2];
    else
        cc = [0.2,0.6,0.7];
    end
    if (nx1==1&&nx2==xl)
        plot([posx(n1)+xbound,posx(n2)],[posy(n1),posy(n2)],'linewidth',max(0.5,20*abs(Vu(i))),'color',cc);
    elseif nx2==1&&nx1==xl
        plot([posx(n2)+xbound,posx(n1)],[posy(n2),posy(n1)],'linewidth',max(0.5,20*abs(Vu(i))),'color',cc);
    else
        plot([posx(n1),posx(n2)],[posy(n1),posy(n2)],'linewidth',max(0.5,20*abs(Vu(i))),'color',cc);
    end
end
%for a = 1:2 % subtract the translational and rotational field
%    displ(nid) = displ(nid) - (disp'*V(:,a))*V(nid,a);
%end
[VV,dd] = eig(full(Mmat));
Vu = VV(:,ff);
quiver(posx',posy',Vu(1:2:end),Vu(2:2:end),'linewidth',2,'AutoScaleFactor',1,'color',[0.1,0.1,0.1])
%partr = sum(displ.^2)^2/num/dim/sum(displ.^4);
xlim([-.5,13])
ylim([-1,11])
