%[pos,nb] = diffNetworks(256,2,5.0,4.0,1);
num = 256;
dim = 2;
z = 4.25;
r = 1;
upp = 0.6;
low = 0.4;
k0  = .1;
% read the network
[pos,nb] = SimpNetwork(num,dim,z,r);
[Smat,bpos] = StructureMat(pos,nb);
% make the bonds in [low,upp] softer
flg = bpos(:,2)>low&bpos(:,2)<upp;
Smat(flg,:) = sqrt(k0)*Smat(flg,:);
%plotNetwork0(pos,nb,1,flg);
%plot(bpos(:,1),bpos(:,2),'bo');
% compute point response
[disp,pert,fext] = pointResponse([.1,0.5],Smat,pos);
%quiver(pos(:,1),pos(:,2),fext(1:2:end),fext(2:2:end),'linewidth',2,'AutoScaleFactor',1,'color',[0.6,0.,0.8])
%quiver(pos(:,1),pos(:,2),disp(1:2:end),disp(2:2:end),'linewidth',2,'AutoScaleFactor',1/sqrt(k0),'color',[0.1,0.1,0.1])
%plot([0.,1.],[upp,upp],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([0.,1.],[low,low],'--','linewidth',2,'color',[0.85,0.4,0]);
%plot([1.,1.],[0,1],'--','linewidth',2,'color',[0.85,0.4,0]);