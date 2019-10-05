function h=mapfigure(posx,posy,z,xb,x0)
%% plot map for quantities
% posx and posy are column vectors of x positions and y positions of data
% sites, z is the column vector of the plotting values on sites
% xb is the size in horizontal direction, x0 is the horizontal shift
cmin = 0; cmax = 1; % range of z value (one can change based on the need)
if nargin>=4 % periodic boundary, reproduce one copy so that the map looks better
    if nargin==4
        x0 = 0;
    end
    posy1 = [posy;posy];
    posx1 = [posx;posx+xb];
    posx1 = rem(posx1+(xb+x0)/2,2*xb);
    zz    = [z;z];
elseif nargin==3
    posx1 = posx;
    posy1 = posy;
    zz    = z;
else
    error('number of input argument wrong!');
end
h=figure;
tri = delaunay(posx1,posy1); % make a triangulation for plot the map
trisurf(tri,posx1,posy1,zz) % 
view(0, 90);
shading interp
axis equal off
caxis([cmin,cmax]) % colored range of z value
if nargin>=4 % adjust the frame
    ylim([0,xb*sqrt(3)/2])
    xlim([(xb-x0)/2-0.5,xb+(xb-x0)/2-0.5])
end