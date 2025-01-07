close all
clear all
clc

global c_parameter alphaalpha AS academic_test P f miu ni Radius DomInf IfPhi2 Grav_const Rb TRadius HydroSolFile ToughSolFile Software FileTopo Ks expmiu

AS=0;
academic_test=0;

IfPhi2=0;
DomInf=0;
expmiu=1;
alphaalpha=1;

depth=500;
FFF=1;

deformation=1;
gravitational_potential=0;

n=150;
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
    solgp = solgp_(Xu,Yu); % exact solution in matrix fashion
    
    % gamma coefficient
    G=G_(Xu,Yu);
    NU=nu_(Xu,Yu);
    
else
    % % % computing the signed distance function from a level-set function
    Phi1=phi1__(Xu,Yu);
    %Phi1=phi1Topo(Xu,Yu,FileTopo);
    if DomInf
        Phi1=-1*ones(size(Phi1));
    end
    if IfPhi2
        Phi2=phi2__(Xu,Yu);
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
    
    % gamma coefficient
    G=G__(Xu,Yu);
    NU=nu__(Xu,Yu);
    Kbulk_s=Kbulk_s__(Xu,Yu);
end

G=miu*ones(size(G));
NU=0.25*ones(size(NU));

% Phi1=Yu-depth;

%%% dimension control: if the arguments have length 1, then transorm it in constant matrix with dimension nn*mm
[Phi1,Phi2,G,NU,Xt1inv,Yt1inv]=controlDim(nn,mm,Phi1,Phi2,G,NU,Xt1inv,Yt1inv);

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

if academic_test
    solP = solp_(Xu,Yu); % exact solution in matrix fashion
    solT = solT_(Xu,Yu); % exact solution in matrix fashion
    Alpha = alpha_(Xu,Yu);
    Beta = beta_(Xu,Yu);
else
    solP=zeros(nn,mm);
    solT=zeros(nn,mm);
    dRho=zeros(nn,mm);
    % % %     if strcmp(Software,'hydro')
    % % %         [solP solT dRho]=Read_hydro_sol(HydroSolFile,Xu,Yu,iter,steady);
    % % %     elseif strcmp(Software,'tough')
    % % %         [solP solT dRho]=Read_TOUGH_sol(ToughSolFile,Xu,Yu,iter,steady);
    % % % %         [solP solT dRho]=Read_TOUGH_sol(ToughSolFile,Xu,Yu,steady,steady); %if you want to consider only the magmatic chamber
    % % %     elseif strcmp(Software,'mufits')
    % % %         [solP solT dRho]=Read_MUFITS_sol(MufitsSolFile,Xu,Yu,iter,steady);
    % % %     end
    [Rx0,Ry0]=RxRy_time(0);
    Phi20=sqrt(((Xu)/Rx0).^2+((Yu+f)/Ry0).^2)-1;
    %     dRho(Phi2<0 & Phi20>0)=-RhoMagma;
    %     dRho(Phi2<0)=-RhoMagma;
    Alpha = alpha__(Xu,Yu);
    Beta = beta__(Xu,Yu);
end

%%% dimension control: if the arguments have length 1, then transorm it in constant matrix with dimension nn*mm
[Alpha,Beta]=controlDim(nn,mm,Alpha,Beta);

Kbulk = 2*G.*(1+NU)/3./(1-2*NU);
Tterm = Kbulk.*Alpha.*solT;
% Tterm = 5e9.*Alpha.*solT;
% Tterm = Kbulk.*Alpha.*solT;
% Pterm = Beta.*(1-Kbulk./Ks).*solP;
% Pterm = Beta.*(1-5./30).*solP;
Pterm = Beta.*solP;
% Pterm = (1-Kbulk./Kbulk_s).*solP;

% Pterm=0*Pterm;
% Tterm=0*Tterm;

I=2:nn-1;
J=2:mm-1;
dPx = zeros(nn,mm);
dPy = zeros(nn,mm);
dTx = zeros(nn,mm);
dTy = zeros(nn,mm);
dPx(I,J) = (Pterm(I,J+1)-Pterm(I,J-1))/2./dx;
dPy(I,J) = (Pterm(I+1,J)-Pterm(I-1,J))/2./dy;
dTx(I,J) = (Tterm(I,J+1)-Tterm(I,J-1))/2./dx;
dTy(I,J) = (Tterm(I+1,J)-Tterm(I-1,J))/2./dy;

% source term F
% [Fu,Fv]=GradSource(Xu,Yu,5000,1000,10e6);
F.u = academic_test*fu_(Xu,Yu); %-Xt1inv.*(dPx+dTx); %+Fu;
F.v = academic_test*fv_(Xu,Yu); %-Yt1inv.*(dPy+dTy); %+Fv;
%dRho=drho__(Xu,Yu);
F.gp = academic_test*fgp_(Xu,Yu)+(~academic_test)*4*pi*Grav_const*dRho;

%%% setting the Grid (Inside/Ghost) structure
[Ghost]=setGrid_P_and_T(Ghost,Pterm,Tterm);

if academic_test
    %%% setting the BC
    [F] = setBC(F,Ghost,BC);
    [F] = setBC_ell(F,Xu,Yu);
else
    if ~DomInf
        F.u([Ghost.Phi1.index])=0;
        F.v([Ghost.Phi1.index])=0;
    end
    
    if IfPhi2
        F.u([Ghost.Phi2.index])=-P*Nx2;
        F.v([Ghost.Phi2.index])=-P*Ny2;
    end
    
    F.u([Ghost.Bdy.index])=0;
    F.v([Ghost.Bdy.index])=0;
    F.gp([1,nn],2:mm-1)=0;
    F.gp(2:nn-1,[1,mm])=0;
end

if ~DomInf
    F.u([Ghost.Phi1.index])=F.u([Ghost.Phi1.index]);%+*([Ghost.Phi1.Tterm]+[Ghost.Phi1.Pterm]).*Nx1;
    F.v([Ghost.Phi1.index])=F.v([Ghost.Phi1.index]);%+*([Ghost.Phi1.Tterm]+[Ghost.Phi1.Pterm]).*Ny1;
end
if IfPhi2
    F.u([Ghost.Phi2.index])=F.u([Ghost.Phi2.index])+([Ghost.Phi2.Tterm]+[Ghost.Phi2.Pterm]).*Nx2;
    F.v([Ghost.Phi2.index])=F.v([Ghost.Phi2.index])+([Ghost.Phi2.Tterm]+[Ghost.Phi2.Pterm]).*Ny2;
end

F.u=zeros(nn,mm);
F.v=zeros(nn,mm);
x_source=0;
y_source=-1*depth;
RR=sqrt((Xu-x_source).^2+(Yu-y_source).^2);
epsilon=1*dx*c_parameter;
Delta_dirac=-exp(-RR.^2/2/epsilon^2)/epsilon^2/(2*pi);
I=2:nn-1;
J=2:mm-1;
F.v(I,J) = Xt1inv(I,J).*(Delta_dirac(I,J+1)-Delta_dirac(I,J-1))/2./dx;
% F.v(I,J) = Yt1inv(I,J).*(Delta_dirac(I+1,J)-Delta_dirac(I-1,J))/2./dy;

% yu=yt__(y);
% [~,row_source]=min(abs(yu+depth));
% F.u(row_source,n/2+1-1)=-FFF;
% F.u(row_source,n/2+1+1)=FFF;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if AS
    if deformation
        [u,v,M,RHS,U]=LinearSystemSparse2hetAS(dx,dy,R,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
    end
    if gravitational_potential
        [gp,M,RHS,GP]=LinearSystemSparseEllipticAS(dx,dy,R,Inside,Ghost,F,Xt1inv,Yt1inv);
        dG=zeros(nn,mm);
        %         Ins_all=find(X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end));
        %         dG(Ins_all)=-(gp(Ins_all+1)-gp(Ins_all-1))/2/dy;
        Ins_all_y=find(Y>Y(1) & Y<Y(end));
        dG(Ins_all_y)=-(gp(Ins_all_y+1)-gp(Ins_all_y-1))/2/dy;
        
        dG=dG.*Yt1inv;
    end
else
    if deformation
        [u,v,M,RHS,U]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
    end
end
% [u_LS,v_LS,M,RHS,U,Pinv,B]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv);
% [u,v,res_u,res_v]=PoissonSolver(x,y,dx,dy,Phi1,Phi2,F,G,NU,BC,Xt1inv,Yt1inv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

Xsup=[Ghost.Phi1.xc];
Ysup=[Ghost.Phi1.yc];
Xusup=xt__(Xsup);
Yusup=yt__(Ysup);

usup=zeros(1,length(Ghost.Phi1));
vsup=zeros(1,length(Ghost.Phi1));
for ii=1:length(Ghost.Phi1)
    usup(ii)=u(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    vsup(ii)=v(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
end
usup=usup';
vsup=vsup';

return;

if deformation
    if ~academic_test
        
        if ~DomInf
            Xsup=[Ghost.Phi1.xc];
            Ysup=[Ghost.Phi1.yc];
            Xusup=xt__(Xsup);
            Yusup=yt__(Ysup);
            [Ux1,Uz1]=MogiSup(Xusup',Yusup',Radius,f,P,miu,ni);
            
            usup=zeros(1,length(Ghost.Phi1));
            vsup=zeros(1,length(Ghost.Phi1));
            for ii=1:length(Ghost.Phi1)
                usup(ii)=u(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
                vsup(ii)=v(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
            end
            usup=usup';
            vsup=vsup';
            
            uchamber=zeros(1,length(Ghost.Phi2));
            vchamber=zeros(1,length(Ghost.Phi2));
            for ii=1:length(Ghost.Phi2)
                uchamber(ii)=u(Ghost.Phi2(ii).stencil(:))'*Ghost.Phi2(ii).coeffsD(:);
                vchamber(ii)=v(Ghost.Phi2(ii).stencil(:))'*Ghost.Phi2(ii).coeffsD(:);
            end
            uchamber=uchamber';
            vchamber=vchamber';
            
            Interp_u=scatteredInterpolant(Xu(Inside.all),Yu(Inside.all),u(Inside.all));
            Interp_v=scatteredInterpolant(Xu(Inside.all),Yu(Inside.all),v(Inside.all));
            u_axis=Interp_u(0,0);
            v_axis=Interp_v(0,0);
        else
            prendi=find(Yu==0);
            usup=u(prendi);
            vsup=v(prendi);
            Xusup=Xu(prendi);
            Yusup=Yu(prendi);
            %[Ux1,Uz1]=DomInfSol(Xusup,Yusup,Radius,f,P,miu,ni);
            [Ux1,Uz1]=Sol_An_Temp(Xusup,Yusup,Radius,Rb,TRadius,miu);
        end
        
        %         figure(1)
        %         subplot(1,2,1)
        %         plot(Xusup,[usup Ux1])
        %         legend('Numeric','Analytic')
        %         xlim([0 10000])
        %         subplot(1,2,2)
        %         plot(Xusup,[vsup -Uz1])
        %         legend('Numeric','Analytic')
        %         xlim([0 10000])
        %         return;
    else
        
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
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     Xsup=[Ghost.Phi1.xc];
    %     Ysup=[Ghost.Phi1.yc];
    %     Xusup=xt_(Xsup);
    %     Yusup=yt_(Ysup);
    %     usup=zeros(1,length(Ghost.Phi1));
    %     vsup=zeros(1,length(Ghost.Phi1));
    %     for ii=1:length(Ghost.Phi1)
    %         usup(ii)=u(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    %         vsup(ii)=v(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    %         solusup(ii)=solu(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    %         solvsup(ii)=solv(Ghost.Phi1(ii).stencil(:))'*Ghost.Phi1(ii).coeffsD(:);
    %     end
    %
    %     plot(Xusup,[usup' solusup'],'*')
    %     figure(100)
    %     plot(Xusup,[vsup' solvsup'],'*')
    %     drawnow;
    %
    %     p_norm=1;
    %     ERR=err_u(:);
    %     L1Su=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    %     p_norm=inf;
    %     LiSu=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    %     disp(['L^1 error u: ',num2str(L1Su),', L^inf error u: ',num2str(LiSu)])
    %     p_norm=1;
    %     ERR=err_v(:);
    %     L1Sv=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    %     p_norm=inf;
    %     LiSv=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    %     disp(['L^1 error v: ',num2str(L1Sv),', L^inf error v: ',num2str(LiSv)])
    
end
if gravitational_potential
    if ~academic_test
        prendi=find(Yu==0);
        dGsup=dG(prendi);
        Xusup=Xu(prendi);
        Yusup=Yu(prendi);
        GravAn=GravSol(Xusup,Yusup,Radius,f,200);%dG in miuGal and drho in kg/m^3
        
        % % %         traccia=load(FileTopo);
        % % %         dGsupT=griddata(Xu(Inside.all),Yu(Inside.all),dG(Inside.all),traccia(:,1),traccia(:,2));
        
        % % %         figure(4)
        % % %         plot(Xusup,[dGsup*1e8 -GravAn],'*', traccia(:,1),[dGsupT*1e8],'ro') %in miuGal
        %plot(Xusup,[dGsup*1e8],'*', traccia(:,1),[dGsupT*1e8],'ro')
        legend('Numeric','Analytic', 'onTopo')
        xlim([0 10000])
        
        %         dGforGrad=griddata(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),dG(2:end-1,2:end-1),Xusup(2:end-1),vsup);
        return;
    end
    Ins_all=find(X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end));
    PTS=setdiff(1:nn*mm,[1 nn nn*mm-nn+1 nn*mm]);
    solgp=zeros(nn,mm);
    solgp(PTS) = solgp_(Xu(PTS),Yu(PTS)); % exact solution in matrix fashion
    SOL=solgp(:);
    
    err_gp=zeros(nn,mm);
    err_gp(Inside.all)=gp(Inside.all)-solgp(Inside.all);
    
    ERR=[err_gp(:)];
    p_norm=1;
    L1S=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    p_norm=inf;
    LiS=norm(ERR,p_norm)/(length(ERR))^(1/p_norm);
    disp(['L^1 error gp: ',num2str(L1S),', L^inf error gp: ',num2str(LiS)])
    
    gpx=zeros(nn,mm);    gpy=zeros(nn,mm);
    gpx(Ins_all)=(gp(Ins_all+nn)-gp(Ins_all-nn))/2/dx;
    gpy(Ins_all)=(gp(Ins_all+1)-gp(Ins_all-1))/2/dy;
    
    gpx=gpx.*Xt1inv;
    gpy=gpy.*Yt1inv;
    
    solgpx=zeros(nn,mm);    solgpy=zeros(nn,mm);
    solgpx(Ins_all)=solgpx_(Xu(Ins_all),Yu(Ins_all));
    solgpy(Ins_all)=solgpy_(Xu(Ins_all),Yu(Ins_all));
    
    PTS=[Ins_all(:)];
    err_gpx=zeros(nn,mm); err_gpy=zeros(nn,mm);
    err_gpx(PTS)=gpx(PTS)-solgpx(PTS);
    err_gpy(PTS)=gpy(PTS)-solgpy(PTS);
    ERRgrad=[sqrt((gpx(PTS)-solgpx(PTS)).^2+(gpy(PTS)-solgpy(PTS)).^2)];
    p_norm=1;
    L1G=norm(ERRgrad,p_norm)/(length(ERRgrad))^(1/p_norm);
    p_norm=inf;
    LiG=norm(ERRgrad,p_norm)/(length(ERRgrad))^(1/p_norm);
    disp(['L^1 error grad gp: ',num2str(L1G),', L^inf error grad gp: ',num2str(LiG)])
end
