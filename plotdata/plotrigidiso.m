
x0 = 0.85;
r  = 3;
d  = 0.;
xl = 2;
[pos,nb,fl,fr] = rigidiso2(64,2,5.,4.,x0,r,d); % cutNetwork(16,16,2,4.25,1,x0); %pos(:,1) = rem(pos(:,1)+1-x0,1);
%[pos,nb] = diffNetworks(64,2,5.0,4.0,0.5,1);

[Smat,bpos] = StructureMat(pos,nb);
%
[disp,pert,fext] = pointResponse([0.0,0.5],Smat,pos);
id  = find(pert~=0);
nid = setdiff(1:numel(pos),id);
%{
%fl = [18,29,36,40]; %
id = [2*fl-1;2*fl];
pert = zeros(numel(pos),1);
pert(id) = randn(2*length(fl),1);
pert(id(1:2:end)) = pert(id(1:2:end))-mean(pert(id(1:2:end)));
pert(id(2:2:end)) = pert(id(2:2:end))-mean(pert(id(2:2:end)));
pert = pert/norm(pert);
nid = setdiff(1:numel(pos),id);
Qmat = zeros(numel(pos));
for i=id
    Qmat(i,i) = 1;
end
Mmat = Smat'*Smat;
Qmat(:,nid) = -Mmat(:,nid);
fext = Mmat*pert;
disp = pinv(Qmat)*fext;
%}
displ = zeros(numel(pos),1);
xid = mod(nid,2)==1;
yid = mod(nid,2)==0;
dx = mean(disp(nid(xid)));
dy = mean(disp(nid(yid)));
displ(nid) = disp(nid)-(dx)*xid'-(dy)*yid';
%{
    y0 = mean(pos(:,1));%.*abs(displ(2:2:end)))/mean(abs(displ(2:2:end)));
    nm = norm(pos(:,1)-y0);
    dx1 = (pos(:,1)-y0)'*displ(2:2:end)/nm;
    displ(2:2:end) = displ(2:2:end)-dx1*(pos(:,1)-y0)/nm;
%
r0 = mean(pos); % position of the rotation center
vr = zeros(numel(pos),1);
vr(1:2:end) = -(pos(:,2)-r0(2));
vr(2:2:end) = pos(:,1)-r0(1);
ar = vr'*displ/norm(vr)^2;
displ = displ-ar*vr;
%}
displ(id) = 0;
%sites = findIsostatic(fr(ceil(4*fr/size(pos,1))==2)',fr,nb,fliplr(pos));
%if ~isempty(sites)
%    plot(pos(sites,2),pos(sites,1),'o','linewidth',2,'color',[0.3,0.9,0],'MarkerFaceColor',[0.33,1,0],'MarkerSize',10)
%    plot(pos(sites,2)-xl,pos(sites,1),'o','linewidth',2,'color',[0.3,0.9,0],'MarkerFaceColor',[0.33,1,0],'MarkerSize',10)
%end
%% plot network
plot(pos(:,2),pos(:,1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',6)
axis equal off
hold all
plot(pos(:,2)-xl,pos(:,1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',6)
for n = 1:size(nb,1)
    n1 = nb(n,1);
    n2 = nb(n,2);
    dx = pos(n2,1)-pos(n1,1);
    dy = pos(n2,2)-pos(n1,2);
    if n1>64&&n1<=128&&n2>64&&n2<=128%bpos(n,1)>1.0 && bpos(n,2)<1 && bpos(n,2)>0. %ismember(n1,sites)&&ismember(n2,sites)%
        cc = [0.2,0.9,0.1];%[0.1,0.2,0.9];%[0.2,0.8,0.1];
    else
        cc = [0.9,0.1,0.2];
    %else
    %    cc = [0.1,0.2,0.9];
    end
    if abs(dx)<0.5 && abs(dy)<0.5
        %putSpring([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],theta,[120/255 80/255 0],1.5);
        %if (pos(n1,1)-x0)*(pos(n2,1)-x0)<=0 || (pos(n1,1)-1+x0)*(pos(n2,1)-1+x0)<=0 || ...
        %        (pos(n1,2)-x0)*(pos(n2,2)-x0)<0 || (pos(n1,2)-1+x0)*(pos(n2,2)-1+x0)<0
        %    plot([pos(n1,1) pos(n1,1)+dx],[pos(n1,2) pos(n1,2)+dy],'--','color',[0 10/255 240/255],'linewidth',3);
        %else
        plot([pos(n1,2) pos(n1,2)+dy],[pos(n1,1) pos(n1,1)+dx],'-','color',cc,'linewidth',3);
        plot([pos(n1,2)-xl pos(n1,2)-xl+dy],[pos(n1,1) pos(n1,1)+dx],'-','color',cc,'linewidth',3);
        %enc
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
        plot([pos(n1,2) pos(n1,2)+dy],[pos(n1,1) pos(n1,1)+dx],'-','color',cc,'linewidth',3);
        plot([pos(n2,2) pos(n2,2)-dy],[pos(n2,1) pos(n2,1)-dx],'-','color',cc,'linewidth',3);
        plot([pos(n1,2)-xl pos(n1,2)-xl+dy],[pos(n1,1) pos(n1,1)+dx],'-','color',cc,'linewidth',3);
        plot([pos(n2,2)-xl pos(n2,2)-xl-dy],[pos(n2,1) pos(n2,1)-dx],'-','color',cc,'linewidth',3);
    end
end
%% vector field

%plot([1.,1.],[1,2],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([0.,0.],[1,2],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([0.,1.],[1,1],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([0.,1.],[2,2],'--','linewidth',2,'color',[0.85,0.4,0]);
[~,ord] = sort(abs(displ),'descend');
%id = ord(1:2);
%displ(135)=0; displ(136)=0;
%displ(id(ceil(id/2)>64))=displ(id(ceil(id/2)>64))/100;
quiver([pos(:,2);pos(:,2)-xl],[pos(:,1);pos(:,1)],[pert(2:2:end);pert(2:2:end)],[pert(1:2:end);pert(1:2:end)],'linewidth',3,'AutoScaleFactor',1,'color',[0.6,0.,0.8])
%quiver([pos(:,2);pos(:,2)-xl],[pos(:,1);pos(:,1)],[fext(2:2:end);fext(2:2:end)],[fext(1:2:end);fext(1:2:end)],'linewidth',2,'AutoScaleFactor',1,'color',[0.8,0.7,0.])
quiver([pos(:,2);pos(:,2)-xl],[pos(:,1);pos(:,1)],[displ(2:2:end);displ(2:2:end)],[displ(1:2:end);displ(1:2:end)],'linewidth',2,'AutoScale','on','AutoScaleFactor',2,'color',[0.1,0.1,0.1])
xlim([-0.77,1.73])
ylim([-0.1,2.2])