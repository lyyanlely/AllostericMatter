%% load the field data in V, cmin, cmax should be adapted to range of V, Lmin, and Lmax are the range of frame
V = zloc{1}; %log10(sbfc(9,:)); %displ(3:3:end); %magd{1}; %
cmin = 5;
cmax = 10;
Lmax = 7.9; %5.1
Lnum = 80;
ticks = cmin:0.5:cmax;
tlab  = ['10^{0.0}';'10^{0.5}';'10^{1.0}';'10^{1.5}';'10^{2.0}';'10^{2.5}'];

% build the 3D frame, and interpolate the field values to the 3D frame
[Xq,Yq,Zq] = meshgrid(-0.1:0.3:Lmax+0.1,-0.1:0.3:Lmax+0.1,-0.1:0.3:Lmax+0.1);
Vt = griddata(pos(1,:),pos(2,:),pos(3,:),V,Xq,Yq,Zq,'nearest');
Vq = griddata(pos(1,:),pos(2,:),pos(3,:),V,Xq,Yq,Zq);
Vq(isnan(Vq)) = Vt(isnan(Vq));
slice(Xq,Yq,Zq,Vq,4.,[2.6,4.5,6.5],7.6); %3.,.2,[0.1,Lmax-0.1]
shading interp
caxis([cmin,cmax]);
%axis equal
box on
xlim([0,Lmax+0.1])
ylim([1.5,7.5]);
zlim([-0.2,Lmax+0.1]);

%% 3D gif scanning over one of the axis
% XY sections
surf(permute(Xq(:,:,1),[1,2,3]),permute(Yq(:,:,1),[1,2,3]),permute(Vq(:,:,1),[1,2,3]))
% XZ sections
%surf(permute(Xq(1,:,:),[2,3,1]),permute(Zq(1,:,:),[2,3,1]),permute(Vq(1,:,:),[2,3,1]))
% YZ sections
%surf(permute(Yq(:,1,:),[1,3,2]),permute(Zq(:,1,:),[1,3,2]),permute(Vq(:,1,:),[1,3,2]))
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'FontSize',20);
shading flat
view(0, 90);
axis equal off
colormap hot
caxis([cmin,cmax]);
colorbar;
%colorbar('Ticks',ticks,'TickLabels',tlab);
xlim([0,Lmax]);
ylim([0,Lmax]);
xlabel('$x$');
ylabel('$y$');
ax = gca;
ax.Units = 'pixels';
pos1 = ax.Position;
marg = 30;
rect = [-marg, -marg, pos1(3)+3*marg, pos1(4)+2*marg];
f = getframe(gca,rect);
[im,map] = rgb2ind(f.cdata,256,'nodither');
for zz = 1:ceil(Lmax/0.1)+1
    hold off
    % XY
    surf(permute(Xq(:,:,zz),[1,2,3]),permute(Yq(:,:,zz),[1,2,3]),permute(Vq(:,:,zz),[1,2,3]))
    % XZ
    %surf(permute(Xq(zz,:,:),[2,3,1]),permute(Zq(zz,:,:),[2,3,1]),permute(Vq(zz,:,:),[2,3,1]))
    % YZ
    %surf(permute(Yq(:,zz,:),[1,3,2]),permute(Zq(:,zz,:),[1,3,2]),permute(Vq(:,zz,:),[1,3,2]))
    hold on
    set(gca,'FontSize',20);
    shading flat
    view(0, 90);
    axis equal off
    colormap hot
    caxis([cmin,cmax]);
    colorbar;
    %colorbar('Ticks',ticks,'TickLabels',tlab);
    xlim([0,Lmax]);
    ylim([0,Lmax]);
    xlabel('$x$'); % change accordingly
    ylabel('$y$');
    f=getframe(gca,rect);
    im(:,:,1,zz) = rgb2ind(f.cdata,map,'nodither');
end
%% unfolded sections plot
hold all
% six surfaces of the cube
surf(permute(Xq(Lnum-2,:,:),[2,3,1]),permute(Zq(Lnum-2,:,:),[2,3,1]),permute(Vq(Lnum-2,:,:),[2,3,1]))
surf(2*Lmax-permute(Yq(:,Lnum-2,:),[1,3,2]),permute(Zq(:,Lnum-2,:),[1,3,2]),permute(Vq(:,Lnum-2,:),[1,3,2]))
surf(permute(Yq(:,3,:)-Lmax,[1,3,2]),permute(Zq(:,3,:),[1,3,2]),permute(Vq(:,3,:),[1,3,2]))
surf(permute(Xq(:,:,Lnum-2),[1,2,3]),2*Lmax-permute(Yq(:,:,Lnum-2),[1,2,3]),permute(Vq(:,:,Lnum-2),[1,2,3]))
surf(permute(Xq(:,:,3),[1,2,3]),permute(Yq(:,:,3)-Lmax,[1,2,3]),permute(Vq(:,:,3),[1,2,3]))
surf(permute(Xq(3,:,:),[2,3,1]),3*Lmax-permute(Zq(3,:,:),[2,3,1]),permute(Vq(3,:,:),[2,3,1]))
shading flat
set(gca,'FontSize',20);
view(0, 90);
axis equal off
colormap hot
caxis([cmin,cmax]);
colorbar
xlim([-Lmax,2*Lmax]);
ylim([-Lmax,3*Lmax]);
scatter3(pos(1,it),2*Lmax-pos(2,it),[10,10,10,10],100,'x','LineWidth',4,'MarkerEdgeColor',[0.1,0.6,0.9])
scatter3(pos(1,id),pos(2,id)-Lmax,[10,10,10,10],100,'x','LineWidth',4,'MarkerEdgeColor',[0.6,0,0.8])

%% Slice plot
hold all
% choose your slices, the example below shows the first y slice
surf(permute(Xq(:,1,:),[1,3,2]), permute(Yq(:,1,:),[1,3,2]), permute(Zq(:,1,:),[1,3,2]), permute(Vq(:,1,:),[1,3,2]))
shading flat
set(gca,'FontSize',20);
% choose your view point
view(0, 90);
axis equal off
colormap hot
caxis([cmin,cmax]);
colorbar
% location of perturbation
scatter3(pos(1,it),pos(2,it),[10,10,10,10],100,'x','LineWidth',4,'MarkerEdgeColor',[0.1,0.6,0.9])
scatter3(pos(1,id),pos(2,id),[10,10,10,10],100,'x','LineWidth',4,'MarkerEdgeColor',[0.6,0,0.8])

%% 3D vector field scanning (adapt in the way similar to scalar scanning)
ratio = 0.2;
quiver(permute(Xq(:,:,1),[1,2,3]),permute(Yq(:,:,1),[1,2,3]),ratio*permute(Vx(:,:,1),[1,2,3]),ratio*permute(Vy(:,:,1),[1,2,3]),'autoscale','off','linewidth',1.,'color',[0.1,0.1,0.1]);
set(0,'DefaultTextInterpreter', 'latex');
set(gca,'FontSize',20);
axis equal
xlim([-0.1,Lmax+0.1]);
ylim([-0.1,Lmax+0.1]);
xlabel('$x$');
ylabel('$y$');
ax = gca;
ax.Units = 'pixels';
pos1 = ax.Position;
marg = 30;
rect = [-marg, -1.5*marg, pos1(3)+2*marg, pos1(4)+2*marg];
f = getframe(gca,rect);
[im,map] = rgb2ind(f.cdata,256,'nodither');
for zz = 1:size(Xq,3)
    quiver(permute(Xq(:,:,zz),[1,2,3]),permute(Yq(:,:,zz),[1,2,3]),ratio*permute(Vx(:,:,zz),[1,2,3]),ratio*permute(Vy(:,:,zz),[1,2,3]),'autoscale','off','linewidth',1.,'color',[0.1,0.1,0.1]);
    set(gca,'FontSize',20);
    axis equal
    xlim([-0.1,Lmax+0.1]);
    ylim([-0.1,Lmax+0.1]);
    xlabel('$x$');
    ylabel('$y$');
    f=getframe(gca,rect);
    im(:,:,1,zz) = rgb2ind(f.cdata,map,'nodither');
end