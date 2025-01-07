close all
clear all
clc

global AS academic_test P;

AS=1;
academic_test=1;

n=32;
if academic_test
    createMatlabFunctions(n);
end
[tank,BC,n]=data_symb_n(n);

m=n*(tank(2)-tank(1))/(tank(4)-tank(3));

% Computational domain [a,b]x[c,d]
a=tank(1);
b=tank(2);
c=tank(3);
d=tank(4);

dx=(b-a)/m; %spatial step in the x-diricetion
x=a:dx:b;
dy=(d-c)/n; %spatial step in the y-diricetion
y=c:dy:d;

[X,Y] = meshgrid(x,y);
[nn,mm]=size(X);

Xu=xt__(X);
Yu=yt__(Y);
xt1=xt1__(x);
yt1=yt1__(y);
[Xt1,Yt1] = meshgrid(xt1,yt1);
Xt1inv=1./Xt1;
Yt1inv=1./Yt1;
R=Xu;

if academic_test
    % % % computing the signed distance function from a level-set function
    Phi1=phi1_(Xu,Yu);
    Phi2=phi2_(Xu,Yu);
    
    % exact solution
    solu = solu_(Xu,Yu); % exact solution in matrix fashion
    solv = solv_(Xu,Yu); % exact solution in matrix fashion
    
    % gamma coefficient
    G=G_(Xu,Yu);
    NU=nu_(Xu,Yu);
    
    % source term F
    F.u = 1*fu_(Xu,Yu);
    F.v = 1*fv_(Xu,Yu);
else
    % % % computing the signed distance function from a level-set function
    Phi1=phi1__(Xu,Yu);
    Phi2=phi2__(Xu,Yu);
    
    % exact solution
    solu = solu__(Xu,Yu); % exact solution in matrix fashion
    solv = solv__(Xu,Yu); % exact solution in matrix fashion
    
    % gamma coefficient
    G=G__(Xu,Yu);
    NU=nu__(Xu,Yu);
    
    % source term F
    F.u = 0;
    F.v = 0;
end

%%% dimension control: if the arguments have length 1, then transorm it in constant matrix with dimension nn*mm
[Phi1,Phi2,F.u,F.v,G,NU,Xt1inv,Yt1inv]=controlDim(nn,mm,Phi1,Phi2,F.u,F.v,G,NU,Xt1inv,Yt1inv);

% normals
[Nx.Phi1,Ny.Phi1]=normals((1:nn*mm)',Phi1,dx,dy);
[Nx.Phi2,Ny.Phi2]=normals((1:nn*mm)',Phi2,dx,dy);
Nx.Phi1=reshape(Nx.Phi1,nn,mm); Ny.Phi1=reshape(Ny.Phi1,nn,mm);
Nx.Phi2=reshape(Nx.Phi2,nn,mm); Ny.Phi2=reshape(Ny.Phi2,nn,mm);

%%% setting the Grid (Inside/Ghost) structure
[Inside,Ghost,Mask]=setGrid(Phi1,Phi2,x,y,dx,dy,G,NU,Nx,Ny,Xt1inv,Yt1inv);

if academic_test
    %%% setting the BC
    [F] = setBC(F,Ghost,BC);
else
    F.u([Ghost.Phi1.index])=0;
    F.v([Ghost.Phi1.index])=0;
    Nx2=[Ghost.Phi2.nx].*[Ghost.Phi2.xt1inv];
    Ny2=[Ghost.Phi2.ny].*[Ghost.Phi2.yt1inv];
    Modl=sqrt(Nx2.^2+Ny2.^2);
    F.u([Ghost.Phi2.index])=-P*Nx2./Modl;
    F.v([Ghost.Phi2.index])=-P*Ny2./Modl;
    F.u([Ghost.Bdy.index])=0;
    F.v([Ghost.Bdy.index])=0;
end

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AS
    [u,v,M,RHS,U,Pinv,B]=LinearSystemSparse2hetAS(dx,dy,R,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
else
    [u,v,M,RHS,U,Pinv,B]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
end
% [u_LS,v_LS,M,RHS,U,Pinv,B]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
% [u,v,res_u,res_v]=PoissonSolver(x,y,dx,dy,Phi1,Phi2,F,G,NU,BC,Xt1inv,Yt1inv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

if ~academic_test
    return;
end

PTS=[Inside.all(:);[Ghost.Phi1.index]';[Ghost.Phi2.index]';[Ghost.Bdy.index]'];
solu=zeros(nn,mm);
solv=zeros(nn,mm);
solu(PTS) = solu_(Xu(PTS),Yu(PTS)); % exact solution in matrix fashion
solv(PTS) = solv_(Xu(PTS),Yu(PTS)); % exact solution in matrix fashion
SOL=[solu(:);solv(:)];

err_u=zeros(nn,mm);
err_v=zeros(nn,mm);
err_u(Inside.all)=u(Inside.all)-solu(Inside.all);
err_v(Inside.all)=v(Inside.all)-solv(Inside.all);

ERR=[err_u(:);err_v(:)];
p_norm=1;
L1S=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
p_norm=inf;
LiS=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
disp(['L^1 error: ',num2str(L1S),', L^inf error: ',num2str(LiS)])

ux=zeros(nn,mm);    uy=zeros(nn,mm);
vx=zeros(nn,mm);    vy=zeros(nn,mm);
ux(Inside.all)=(u(Inside.all+nn)-u(Inside.all-nn))/2/dx;
uy(Inside.all)=(u(Inside.all+1)-u(Inside.all-1))/2/dy;
vx(Inside.all)=(v(Inside.all+nn)-v(Inside.all-nn))/2/dx;
vy(Inside.all)=(v(Inside.all+1)-v(Inside.all-1))/2/dy;

ux=ux.*Xt1inv;
uy=uy.*Yt1inv;
vx=vx.*Xt1inv;
vy=vy.*Yt1inv;

solux=zeros(nn,mm);    soluy=zeros(nn,mm);
solux(Inside.all)=solux_(Xu(Inside.all),Yu(Inside.all));
soluy(Inside.all)=soluy_(Xu(Inside.all),Yu(Inside.all));
solvx=zeros(nn,mm);    solvy=zeros(nn,mm);
solvx(Inside.all)=solvx_(Xu(Inside.all),Yu(Inside.all));
solvy(Inside.all)=solvy_(Xu(Inside.all),Yu(Inside.all));

PTS=[Inside.all(:)];
err_ux=zeros(nn,mm); err_uy=zeros(nn,mm);
err_ux(PTS)=ux(PTS)-solux(PTS);
err_uy(PTS)=uy(PTS)-soluy(PTS);
err_vx=zeros(nn,mm); err_vy=zeros(nn,mm);
err_vx(PTS)=vx(PTS)-solvx(PTS);
err_vy(PTS)=vy(PTS)-solvy(PTS);
ERRgrad=[sqrt((ux(PTS)-solux(PTS)).^2+(uy(PTS)-soluy(PTS)).^2);sqrt((vx(PTS)-solvx(PTS)).^2+(vy(PTS)-solvy(PTS)).^2)];
p_norm=1;
L1G=norm(ERRgrad,p_norm)/(length(ERRgrad))^(1/p_norm);
p_norm=inf;
LiG=norm(ERRgrad,p_norm)/(length(ERRgrad))^(1/p_norm);
disp(['L^1 error grad: ',num2str(L1G),', L^inf error grad: ',num2str(LiG)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xsup=[Ghost.Phi1.xc];
Ysup=[Ghost.Phi1.yc];
Xusup=xt_(Xsup);
Yusup=yt_(Ysup);
usup=zeros(1,length(Ghost.Phi1));
vsup=zeros(1,length(Ghost.Phi1));
for ii=1:length(Ghost.Phi1)
    usup(ii)=u(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    vsup(ii)=v(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    solusup(ii)=solu(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    solvsup(ii)=solv(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);    
end

plot(Xusup,[usup' solusup'],'*')
figure(100)
plot(Xusup,[vsup' solvsup'],'*')
drawnow;

p_norm=1;
ERR=err_u(:);
L1Su=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
p_norm=inf;
LiSu=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
disp(['L^1 error u: ',num2str(L1Su),', L^inf error u: ',num2str(LiSu)])
p_norm=1;
ERR=err_v(:);
L1Sv=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
p_norm=inf;
LiSv=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
disp(['L^1 error v: ',num2str(L1Sv),', L^inf error v: ',num2str(LiSv)])
