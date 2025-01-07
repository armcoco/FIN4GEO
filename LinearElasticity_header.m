% close all
% clear all
% clc

global AS academic_test DomInf IfPhi2

AS=1;
if ~exist('academic_test','var') || isempty(academic_test)
    academic_test=0;
end
deformation=1;
gravitational_potential=1;
self_potential=1;

% n=416;
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
RR=struct();
RR.centre=Xu;
I=2:nn-1;
J=2:mm-1;
RR.left=nan(nn,mm);
RR.right=nan(nn,mm);
RR.left(I,J)=xt__(0.5*(X(I,J)+X(I,J-1)));
RR.right(I,J)=xt__(0.5*(X(I,J)+X(I,J+1)));

if academic_test
    % % % computing the signed distance function from a level-set function
    Phi1=phi1_(Xu,Yu);
    Phi2=phi2_(Xu,Yu);
    
    % exact solution
    solu = solu_(Xu,Yu); % exact solution in matrix fashion
    solv = solv_(Xu,Yu); % exact solution in matrix fashion
    solgp = solgp_(Xu,Yu); % exact solution in matrix fashion
    solsp = solsp_(Xu,Yu); % exact solution in matrix fashion
    
    % grav coefficients
    G=G_(Xu,Yu);
    NU=nu_(Xu,Yu);
    % SP coefficients
    Sig=sig_(Xu,Yu);
    L_SP=L_SP_(Xu,Yu);
    
else
    % % % computing the signed distance function from a level-set function
    Phi1=phi1__(Xu,Yu);
    %Phi1=phi1Topo(Xu,Yu,FileTopo);
    if DomInf
        Phi1=-1*ones(size(Phi1));
    end
    if IfPhi2
        Phi2=phi2__(Xu,Yu);
        Phi2 = max(max(Xu-4500,-6500.3-Yu),Yu+6500); %for Gottsmann 23 Jan 2018
    else
        Phi2=1;
    end
    % exact solution
    solu = solu__(Xu,Yu); % exact solution in matrix fashion
    solv = solv__(Xu,Yu); % exact solution in matrix fashion
    %   solP = solp__(Xu,Yu); % exact solution in matrix fashion
    %   solT = solT__(Xu,Yu); % exact solution in matrix fashion
    %   dRho=drho__(Xu,Yu);
    solgp = solgp__(Xu,Yu); % exact solution in matrix fashion
    solsp = solsp__(Xu,Yu); % exact solution in matrix fashion
    
    % grav coefficients
    G=G__(Xu,Yu);
    NU=nu__(Xu,Yu);
%     G=6.6e9*ones(size(G)); %for Gottsmann 23 Jan 2018
%     NU=0.25*ones(size(NU)); %for Gottsmann 23 Jan 2018
    Kbulk_s=Kbulk_s__(Xu,Yu);
    % SP coefficients
    Sig=sig__(Xu,Yu);
    L_SP=L_SP__(Xu,Yu);
end

%%% dimension control: if the arguments have length 1, then transorm it in constant matrix with dimension nn*mm
[Phi1,Phi2,G,NU,Sig,L_SP,Xt1inv,Yt1inv]=controlDim(nn,mm,Phi1,Phi2,G,NU,Sig,L_SP,Xt1inv,Yt1inv);

% normals
[Nx.Phi1,Ny.Phi1]=normals((1:nn*mm)',Phi1,dx,dy);
[Nx.Phi2,Ny.Phi2]=normals((1:nn*mm)',Phi2,dx,dy);
Nx.Phi1=reshape(Nx.Phi1,nn,mm); Ny.Phi1=reshape(Ny.Phi1,nn,mm);
Nx.Phi2=reshape(Nx.Phi2,nn,mm); Ny.Phi2=reshape(Ny.Phi2,nn,mm);

%%% setting the Grid (Inside/Ghost) structure
[Inside,Ghost,Mask]=setGrid(Phi1,Phi2,x,y,dx,dy,G,NU,Nx,Ny,Xt1inv,Yt1inv);

if ~DomInf
Nx1=[Ghost.Phi1.nx].*[Ghost.Phi1.xt1inv];
Ny1=[Ghost.Phi1.ny].*[Ghost.Phi1.yt1inv];
Modl1=sqrt(Nx1.^2+Ny1.^2);
Nx1=Nx1./Modl1;
Ny1=Ny1./Modl1;
end

if IfPhi2
Nx2=[Ghost.Phi2.nx].*[Ghost.Phi2.xt1inv];
Ny2=[Ghost.Phi2.ny].*[Ghost.Phi2.yt1inv];
Modl2=sqrt(Nx2.^2+Ny2.^2);
Nx2=Nx2./Modl2;
Ny2=Ny2./Modl2;
end

