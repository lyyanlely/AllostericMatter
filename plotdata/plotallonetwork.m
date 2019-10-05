function plotallonetwork(config,nc,id,it,posx,posy,bondnb,displ,pert,targ)
% config is the matrix stores configurations, each row stores one
% configuration in zeros and ones, each column labels one specific bond
% position, one for being occupied by spring, and zero for being vacant.
% nc is the index of the configuration to plot
% id is the label of sites in the stimulus
% it is the label of sites in the target
% posx is a row vector of x component of site positions
% posy is a row vector of y component of site positions
% bondnb is a nb by 2 matrix, for each one of nb bonds, the labels of the 2
% sites are connected.
% displ is a column vector of the response field in
% (1x,1y,2x,2y,3x,3y,...,Nx,Ny). 
% pert is a column vector of the perturbation field
% targ is a column vector of the target field
xl = 12;
yl = 12;
sf = 0.4; % scale factor of vectors at perturbation and target
sfr= 1.0; % scale factor of vectors of response field
stat = config(nc,:); 
occ = find(stat>0);
nsp = sum(stat);
plot(posx,posy,'o','linewidth',2,'color',[0.9,0.3,0],'MarkerFaceColor',[1,0.33,0],'MarkerSize',3)
axis equal off
hold all
plot(posx(posx<1)+xbound,posy(posx<1),'o','linewidth',2,'color',[0.9,0.3,0],'MarkerSize',3)
for i = 1:nsp
    n = occ(i);
    n1 = bondnb(n,1);
    n2 = bondnb(n,2);
    ny1 = ceil(n1/xl);
    nx1 = n1 - (ny1-1)*xl;
    ny2 = ceil(n2/xl);
    nx2 = n2 - (ny2-1)*xl;
    cc = [0.9,0.1,0.2];  % bond color
    if (nx1==1&&nx2==xl)
        plot([posx(n1)+xbound,posx(n2)],[posy(n1),posy(n2)],'linewidth',1,'color',cc);
    elseif nx2==1&&nx1==xl
        plot([posx(n2)+xbound,posx(n1)],[posy(n2),posy(n1)],'linewidth',1,'color',cc);
    else
        plot([posx(n1),posx(n2)],[posy(n1),posy(n2)],'linewidth',1,'color',cc);
    end
end

quiver(posx(id),posy(id),pert(2*id-1)',pert(2*id)','linewidth',3,'AutoScaleFactor',sf,'color',[0.6,0.,0.8])

dispx = displ(1:2:end);
dispy = displ(2:2:end);
num   = length(posx);
tid   = setdiff((1:num),id);
quiver(posx(it),posy(it),targ(2*it-1)'+mean((displ(2*it-1)-targ(2*it-1))),targ(2*it)'+mean((displ(2*it)-targ(2*it))),'linewidth',3,'AutoScaleFactor',sf,'color',[0.1,0.6,0.9]) %
quiver(posx(tid)',posy(tid)',dispx(tid),dispy(tid),'linewidth',2,'AutoScaleFactor',sfr,'color',[0.1,0.1,0.1])
% 
xlim([-.5,xl+1])
ylim([-1,sqrt(3)*yl/2+sqrt(3)])