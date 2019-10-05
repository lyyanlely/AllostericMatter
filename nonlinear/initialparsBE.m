function pars = initialparsBE(pos,nb,flg,L,kw,a,b)
% initialize the parameters used in Conjugate Gradient
nnb  = size(nb,1);

pars.x0 = reshape(pos',[],1);  % initial position
pars.pos = pos;                % position
pars.L  = L;                   % system size
pars.nb = nb;                  % indices of particle neighbors
pars.xb = 1.;                  % system dimension in x direction 
pars.yb = pars.xb*sqrt(3)/2;   % system dimension in y
pars.l0 = pars.xb/L;           % lattice constant l0
pars.kw = kw;                  % strength of weak springs
pars.a  = a*pars.l0;           % mismatches of weak springs
pars.b  = b;                   % external field

flag = zeros(nnb,1);
for i=1:nnb
    n1 = nb(i,1);
    n2 = nb(i,2);
    if flg(i)
        dx = pos(n2,1)-pos(n1,1);
        dy = pos(n2,2)-pos(n1,2);
        if dx>pars.xb/2
            dx = dx-pars.xb;
        elseif dx<-pars.xb/2
            dx = dx+pars.xb;
        end
        if dy>pars.yb/2
            dy = dy-pars.yb;
        elseif dy<-pars.yb/2
            dy = dy+pars.yb;
        end
        if dx*dy>eps
            flag(i) = 1;
        elseif dx*dy<-eps
            flag(i) = -1;
        else
            flag(i) = 2;
        end
    end
end
pars.flag = flag;             % indicates the directions of weak springs 1 for same as external field, -1 for contradict with external field, 2 for perpendicular, 0 for not weak