
%% prepare
r0 = 1;   % lattice constant typical size of the spring
ng = 10;  % size of grid
kx = -pi/xbound:(2*pi)/ng/xbound:pi/xbound; %/xbound;
ky = -pi/ybound:(2*pi)/ng/ybound:pi/ybound; %/ybound;
[kmx,kmy] = meshgrid(kx,ky);
ix = length(kx);
iy = length(ky);

Sq = cell(size(kmx));
U  = cell(size(kmx));
S  = cell(size(kmx));
V  = cell(size(kmx));
dd = cell(size(kmx));
for ii = 1:ix
    for jj = 1:iy
        Sq{ii,jj} = Smatrix;
        id = mod(bondnb(n,1),xl)==0 & mod(bondnb(n,2),xl)==1;
        Sq{ii,jj}(id,:) = Smatrix(id,:)*exp(1j*kx(ii)*xbound);
        %Sq{ii,jj}(:,1:2:end) = Smatrix(:,1:2:end).*exp(r0/2*1j*(kx(ii)*Smatrix(:,1:2:end)+ky(jj)*Smatrix(:,2:2:end)));
        %Sq{ii,jj}(:,2:2:end) = Smatrix(:,2:2:end).*exp(r0/2*1j*(kx(ii)*Smatrix(:,1:2:end)+ky(jj)*Smatrix(:,2:2:end)));
    end
end
%% dispersion relation
%% topological order
stat = config(end,:,9);
detS = zeros(size(kmx));
Sinv = cell(size(kmx));
omegaq = zeros(size(kmx));
for ii = 1:ix
    for jj = 1:iy
        [U{ii,jj},S{ii,jj},V{ii,jj}] = svd(full(Sq{ii,jj}(stat>0,:)));
        dd{ii,jj} = diag(S{ii,jj});
        detS(ii,jj) = (det(U{ii,jj})*det(V{ii,jj})*exp(sum(log(dd{ii,jj}(1:end-2)))));
        ddd = length(dd{ii,jj});
        idx = abs(dd{ii,jj})>1e-14;
        Sinv{ii,jj} = V{ii,jj}*sparse(find(idx),find(idx),dd{ii,jj}(idx).^-1,size(S{ii,jj},2),size(S{ii,jj},1))*U{ii,jj}';
        %if ii~=ceil(ix/2) || jj~=ceil(iy/2)
            omegaq(ii,jj) = min(abs(dd{ii,jj}(idx)));
        %end
    end
end

ni = zeros(ng,ng);
for ii = 1:ix-1
    for jj = 1:iy-1
        ni(ii,jj) = sum(diag(Sinv{ii,jj}*Sq{ii+1,jj}(stat>0,:)))...
                    +sum(diag(Sinv{ii+1,jj}*Sq{ii+1,jj+1}(stat>0,:)))...
                    +sum(diag(Sinv{ii+1,jj+1}*Sq{ii,jj+1}(stat>0,:)))...
                    +sum(diag(Sinv{ii,jj+1}*Sq{ii,jj}(stat>0,:)))-4*ddd+8;
        ni(ii,jj) = ni(ii,jj)/2/xl/yl/2/pi/1i;
    end
end
%% complex angle
phi = zeros(ng,ng);
for ii = 1:ix-1
    for jj = 1:iy-1
        phs = unwrap(angle([detS(ii,jj),detS(ii+1,jj),detS(ii+1,jj+1),detS(ii,jj+1),detS(ii,jj)]));
        phi(ii,jj) = (phs(end)-phs(1))/2/pi;
    end
end
[klx,kly] = meshgrid(kx(1:end-1)+(kx(2)-kx(1))/2,ky(1:end-1)+(ky(2)-ky(1))/2);
%%
phx = zeros(1,ix);
for ii = 1:ix
    for jj = 1:iy-1
        phx(ii) = phx(ii)+sum(diag(Sinv{ii,jj}*Sq{ii,jj+1}(stat>0,:)))-ddd+2;
    end
    phx(ii) = phx(ii)+sum(diag(Sinv{ii,iy}*Sq{ii,1}(stat>0,:)))-ddd+2;
    phx(ii) = phx(ii)/2/xl/yl/2/pi/1i;
end