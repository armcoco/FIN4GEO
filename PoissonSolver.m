function [u,v,res_u,res_v]=PoissonSolver(x,y,dx,dy,Phi1,Phi2,F,G,NU,BC,Xt1inv,Yt1inv)

global solu solv;

NiterDown=2;
NiterUp=1;
NVcycles=1e3;
TOLLres=1e-6;

nn=length(y);
mm=length(x);
n=nn-1;
m=mm-1;

% normals
[Nx.Phi1,Ny.Phi1]=normals((1:nn*mm)',Phi1,dx,dy);
[Nx.Phi2,Ny.Phi2]=normals((1:nn*mm)',Phi2,dx,dy);
Nx.Phi1=reshape(Nx.Phi1,nn,mm); Ny.Phi1=reshape(Ny.Phi1,nn,mm);
Nx.Phi2=reshape(Nx.Phi2,nn,mm); Ny.Phi2=reshape(Ny.Phi2,nn,mm);

%%% setting the Grid (Inside/Ghost) structure
[Inside,Ghost,Mask,MaskI]=setGrid(Phi1,Phi2,x,y,dx,dy,G,NU,Nx,Ny,Xt1inv,Yt1inv);

% % Initial data
u=0*ones(nn,mm);
v=0*ones(nn,mm);
u(Inside.all)=1;
v(Inside.all)=1;
[X,Y]=meshgrid(x,y);
u=sin(X).*sin(Y);
v=0*ones(nn,mm);
% u([Ghost.Bdy.index])=1;
% v([Ghost.Bdy.index])=1;
% u([Ghost.Phi1.index])=1;
% v([Ghost.Phi1.index])=1;
% u([Ghost.Phi2.index])=1;
% v([Ghost.Phi2.index])=1;
norm_res_old=1;

% [u,v]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC);

% % % u(Inside.all)=solu(Inside.all);
% % % v(Inside.all)=solv(Inside.all);
% % % u=solu;
% % % v=solv;
% % % [res_u,res_v] = defect(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC);

for iter=1:NVcycles
    [u,v] = Vcycle(u,v,x,y,dx,dy,NiterDown,NiterUp,n,m,Phi1,Phi2,Inside,Ghost,Mask,MaskI,F,G,NU,BC,Nx,Ny,Xt1inv,Yt1inv);
%     [u,v] = Relaxation(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
%     save u1.txt -ascii u
%     save v1.txt -ascii v
    [res_u,res_v] = defect(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
    max_res=max(abs([res_u(:);res_v(:)]))/(1+max([abs(F.u(:));abs(F.v(:))]));
    p_norm=inf;
    norm_res=norm([res_u(:);res_v(:)],p_norm);
    rho = norm_res/norm_res_old;
%     rhohist(iter)=rho;
%     rho=prod(rhohist(6:end))^(1/length(rhohist(6:end)));    

    if norm_res<TOLLres
%         break;
    end
    
    ErrorArray=abs([u(Inside.all)-solu(Inside.all) ; v(Inside.all)-solv(Inside.all)]);
    Li=max(ErrorArray);
    L1=norm(ErrorArray,1)/length(ErrorArray);
    disp(['iter=',num2str(iter),'; rho=',num2str(rho),'; res=',num2str(norm_res),', L^1 error: ',num2str(L1),', L^inf error: ',num2str(Li)]);
    
    norm_res_old=norm_res;
end
