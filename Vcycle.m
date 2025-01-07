function [u,v] = Vcycle(u,v,x,y,dx,dy,NiterDown,NiterUp,n,m,Phi1,Phi2,Inside,Ghost,Mask,MaskI,F,G,NU,BC,Nx,Ny,Xt1inv,Yt1inv)

salvare=0;

gamma = 2;
coarsestLevel = 32; % Numero di intervalli in cui � diviso il dominio nella coarsest grid
NiterCoarsestLevel=10;
esp=10;             % Numero di iterazioni per espansioneR

[X,Y]=meshgrid(x,y);

nn=n+1;
mm=m+1;

if n==coarsestLevel || m==coarsestLevel
%     for iter = 1:NiterCoarsestLevel
%         [u,v] = Relaxation(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
%     end
    [u,v]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
    return;
end

for iter = 1:NiterDown
    [u,v] = Relaxation(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
end

if salvare
save u2.txt u -ascii -double
save v2.txt v -ascii -double
end

nc = n/2;
mc = m/2;
nnc=nc+1;
mmc=mc+1;
dxc = dx*2;
dyc = dy*2;
xc=coarseE(x);
yc=coarseE(y);
Phi1c = coarseE(Phi1);
Phi2c = coarseE(Phi2);

[res_u,res_v] = defect(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);

if salvare
save resu.txt res_u -ascii -double
save resv.txt res_v -ascii -double
end

ALL=find(X>X(1) & X<X(nn*mm) & Y>Y(1) & Y<Y(nn));
expres1=setdiff(ALL(:),[find(Phi1(:)<0);[Ghost.Phi1.index]']);
expres2=setdiff(find(Phi2(:)<0),[Ghost.Phi2.index]');

for iter=1:esp
    res_u = espansioneR(res_u,expres1,Nx.Phi1,Ny.Phi1,1);
    res_u = espansioneR(res_u,expres2,Nx.Phi2,Ny.Phi2,-1);
    res_v = espansioneR(res_v,expres1,Nx.Phi1,Ny.Phi1,1);
    res_v = espansioneR(res_v,expres2,Nx.Phi2,Ny.Phi2,-1);
end

if salvare
save resu_esp.txt res_u -ascii -double
save resv_esp.txt res_v -ascii -double
end

resc_u=zeros(nnc,mmc);
resc_v=zeros(nnc,mmc);

% resc_u = coarseFW(resc_u,res_u,Phi1<0);
% resc_v = coarseFW(resc_v,res_v,Phi1<0);
% resc_u = coarseFW(resc_u,res_u,Phi1>=0);
% resc_v = coarseFW(resc_v,res_v,Phi1>=0);

% resc_u = coarseFW2(resc_u,res_u,Phi1<0,Nx.Phi1,Ny.Phi1);
% resc_v = coarseFW2(resc_v,res_v,Phi1<0,Nx.Phi1,Ny.Phi1);
% resc_u = coarseFW2(resc_u,res_u,Phi1>=0,-Nx.Phi1,-Ny.Phi1);
% resc_v = coarseFW2(resc_v,res_v,Phi1>=0,-Nx.Phi1,-Ny.Phi1);

resc_u = coarseFW(resc_u,res_u,MaskI);
resc_u = coarseFW(resc_u,res_u,Phi1>=0);
resc_u = coarseFW(resc_u,res_u,Phi2<0);
resc_v = coarseFW(resc_v,res_v,MaskI);
resc_v = coarseFW(resc_v,res_v,Phi1>=0);
resc_v = coarseFW(resc_v,res_v,Phi2<0);

if salvare
save resc_u.txt resc_u -ascii -double
save resc_v.txt resc_v -ascii -double
end

resc.u=resc_u;
resc.v=resc_v;

Gc=coarseE(G);
NUc=coarseE(NU);
Xt1invc=coarseE(Xt1inv);
Yt1invc=coarseE(Yt1inv);
Nxc.Phi1=coarseE(Nx.Phi1);  Nyc.Phi1=coarseE(Ny.Phi1);
Nxc.Phi2=coarseE(Nx.Phi2);  Nyc.Phi2=coarseE(Ny.Phi2);
[Insidec,Ghostc,Maskc,MaskIc]=setGrid(Phi1c,Phi2c,xc,yc,dxc,dyc,Gc,NUc,Nxc,Nyc,Xt1invc,Yt1invc);

e_u = zeros(nnc,mmc);
e_v = zeros(nnc,mmc);
for tt=1:gamma
    [e_u,e_v] = Vcycle(e_u,e_v,xc,yc,dxc,dyc,NiterDown,NiterUp,nc,mc,Phi1c,Phi2c,Insidec,Ghostc,Maskc,MaskIc,resc,Gc,NUc,BC,Nxc,Nyc,Xt1invc,Yt1invc);
end

if salvare
save e_u.txt e_u -ascii -double
save e_v.txt e_v -ascii -double
end

[Xc,Yc]=meshgrid(xc,yc);
ALLc=find(Xc>Xc(1) & Xc<Xc(nnc*mmc) & Yc>Yc(1) & Yc<Yc(nnc));
expres1c=setdiff(ALLc(:),[find(Phi1c(:)<0);[Ghostc.Phi1.index]']);
expres2c=setdiff(find(Phi2c(:)<0),[Ghostc.Phi2.index]');

for iter=1:esp
    e_u = espansioneR(e_u,expres1c,Nxc.Phi1,Nyc.Phi1,1);
    e_u = espansioneR(e_u,expres2c,Nxc.Phi2,Nyc.Phi2,-1);
    e_v = espansioneR(e_v,expres1c,Nxc.Phi1,Nyc.Phi1,1);
    e_v = espansioneR(e_v,expres2c,Nxc.Phi2,Nyc.Phi2,-1);
end

if salvare
save e_u_esp.txt e_u -ascii -double
save e_v_esp.txt e_v -ascii -double
end

e_uf = fineE(e_u);
e_vf = fineE(e_v);


if salvare
save e_uf.txt e_uf -ascii -double
save e_vf.txt e_vf -ascii -double
end

u=u+e_uf;
v=v+e_vf;

if salvare
save u_interp.txt u -ascii -double
save v_interp.txt v -ascii -double
end

for iter = 1:NiterUp
    [u,v] = Relaxation(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
end

if salvare
save u_after.txt u -ascii -double
save v_after.txt v -ascii -double
end
