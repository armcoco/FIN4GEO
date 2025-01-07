% clear all
close all

FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);



steady=1000;
% unrest=10;
% endsim=48;

%     pathSol='Faults_Alia/Fault_inj_massflux_per_unit_area_constant/unrest_x1/';
%     pathSol='Faults_3Km/ScenarioB/unrest_x1/';
%     pathSol='Alia_flow.out_original/HighInj/';
%     pathSol='Faults_Alia_same_grid_as_Alia/unrest/';
%     pathSol='../MultigridAxiDefoGrav 18nov14'/Alia_flow.inp_original/HighInj_withoutINCON_on_flow.inp/INCON_from_steadystate for Low_Fumar_1e-14';
%     pathSol='Faults/ScenarioC/unrest_x1/';
%     pathSol='Faults_perm_xi_not_lower/ScenarioB/unrest_x1/';
%     pathSol='Faults_perm_xi_not_lower_eta_1213/ScenarioB/unrest_x1/';
%     pathSol='Faults_larger_damage_zone/ScenarioB/unrest_x1/';
%     pathSol='INCON_from_steadystate for Low_Fumar_1e-14/';
%     pathSol='INCON_from_SteadyState_Fumar_1e-13/';
%     pathSol='HighInj/';
%     pathSol='Alia_HighInj/';
%     pathSol='Faults_paper/ScenarioC/unrest_x1/';
    pathSol='Faults_paper_new/ScenarioB/unrest_x1/';
    pathSol='Faults_paper_new_extended_release_200m/unrest_x1p15_molar_ratio_0p17_Inj2000m_hybrid_2/';

Sol=load([pathSol,'Sol_',num2str(steady)]);
X=Sol(:,1);
Z=Sol(:,3);
T0=Sol(:,4);
P0=Sol(:,5);
SG0=Sol(:,6);
a=find(Z==Z(1));
l=length(a);
m=length(Z)/l;
TT0=reshape(T0,l,m);
PP0=reshape(P0,l,m);
SSGG0=reshape(SG0,l,m);


% time=load('faults_3000_rainfall_only/unrest_and_quiet_rainfall_constant/time');
time=load([pathSol,'time']);
for iter=0:length(time)-1
    time_=time(iter+1)/86400/365.25;
    clf
    Sol=load([pathSol,'Sol_',num2str(iter)]);
    Conn=load([pathSol,'Conn_',num2str(iter)]);
    X=Sol(:,1);
    Z=Sol(:,3);
    T=Sol(:,4);
    
    P=Sol(:,5);
    SG=Sol(:,6);
    
    a=find(Z==Z(1));
    l=length(a);
    m=length(Z)/l;
    XX=reshape(X,l,m);
    ZZ=reshape(Z,l,m);
    TT=reshape(T,l,m);
    PP=reshape(P,l,m);
    SSGG=reshape(SG,l,m);
    
    subplot(3,1,1)
%     I=2:l-1;
%     J=2:m-1;
%     PPx=zeros(size(PP));
%     PPy=zeros(size(PP));
%     PPx(I,J)=(PP(I+1,J)-PP(I-1,J))./(XX(I+1,J)-XX(I-1,J));
%     PPy(I,J)=(PP(I,J+1)-PP(I,J-1))./(ZZ(I,J+1)-ZZ(I,J-1));
    [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
    clabel(c,h);
%     caxis([-0.5 30])
%     caxis([0 5])
    %     xlim([0 1000])
    colorbar
    title(['pore pressure [MPa],   time: ',num2str(time_)],'FontSize',16)
    
    subplot(3,1,2)
    [c,h]=contourf(XX,ZZ,TT-TT0);
    clabel(c,h);
%     caxis([0 360])
    caxis([-15 100])
    %     xlim([0 1000])
    colorbar
    title('temperature [C]','FontSize',16)
    hold on
%     plottaGradT;
    
    subplot(3,1,3)
    meshsize;
    Nx=l;
    Nz=m;
    FLOH_u=-reshape(Conn(1:(l-1)*m,1),l-1,m)';
    FLOH_v=-reshape(Conn(((l-1)*m+1):end,1),l,m-1)';
    FLOF_u=-reshape(Conn(1:(l-1)*m,2),l-1,m)';
    FLOF_v=-reshape(Conn(((l-1)*m+1):end,2),l,m-1)';
    
    FLOH_u=FLOH_u./Zsize/2/pi./X_u;
    FLOH_v=FLOH_v./Xsize/2/pi./X_v;
    FLOF_u=FLOF_u./Zsize/2/pi./X_u;
    FLOF_v=FLOF_v./Xsize/2/pi./X_v;
    
    Ii=2:m-1;
    Ji=2:l-1;
    FLOWH_u=zeros(m,l);
    FLOWH_v=zeros(m,l);
    FLOWH_u(:,Ji)=(FLOH_u(:,Ji-1)+FLOH_u(:,Ji))/2;
    FLOWH_u(:,[1,l])=FLOH_u(:,[1,Nx-1]);
    FLOWH_v(Ii,:)=(FLOH_v(Ii-1,:)+FLOH_v(Ii,:))/2;
    FLOWH_v([1,m],:)=FLOH_v([1,m-1],:);
    FLOWF_u=zeros(m,l);
    FLOWF_v=zeros(m,l);
    FLOWF_u(:,Ji)=(FLOH_u(:,Ji-1)+FLOH_u(:,Ji))/2;
    FLOWF_u(:,[1,l])=FLOH_u(:,[1,Nx-1]);
    FLOWF_v(Ii,:)=(FLOH_v(Ii-1,:)+FLOH_v(Ii,:))/2;
    FLOWF_v([1,m],:)=FLOH_v([1,m-1],:);
    
    
    CH=log10(sqrt(FLOWH_u.^2+FLOWH_v.^2));
    CF=log10(sqrt(FLOWF_u.^2+FLOWF_v.^2));
    hold on
    %     contourf(XX',ZZ',CH,'lineWidth',1);
%     contourf(XX',ZZ',CF,'lineWidth',1);
    contourf(XX,ZZ,SSGG-SSGG0,'lineWidth',1);
    hh = colorbar;
    ylabel(hh, ['log of magnitude of flow arrays'])
    plot([X(1) X(end) X(end) X(1) X(1)],[Z(1) Z(1) Z(end) Z(end) Z(1)])
    %     axis equal
    
    Xplot=100; % it is the number of points in $x$-direction for which we want to plot the vectors of the flow.
    Zplot=15; % it is the number of points in $z$-direction for which we want to plot the vectors of the flow.
    xx=linspace(Xcenter(1),Xcenter(end),Xplot*2+1);
    xx=xx(2:2:end);
    zz=linspace(Zcenter(1),Zcenter(end),Zplot*2+1);
    zz=zz(2:2:end);
    [XXplot,ZZplot]=meshgrid(xx,zz);
    
    %     F1H=interp2(X_u,Z_u,FLOH_u,XXplot,ZZplot,'cubic');
    %     F2H=interp2(X_v,Z_v,FLOH_v,XXplot,ZZplot,'cubic');
    %     modl=sqrt(F1H.^2+F2H.^2); %./log(sqrt(F1.^2+F2.^2));
    %     lh=quiver3(XXplot,ZZplot,10*ones(size(XXplot)),F1H./modl,F2H./modl,zeros(size(XXplot)));
    
    F1F=interp2(X_u,Z_u,FLOF_u,XXplot,ZZplot,'cubic');
    F2F=interp2(X_v,Z_v,FLOF_v,XXplot,ZZplot,'cubic');
    modl=sqrt(F1F.^2+F2F.^2); %./log(sqrt(F1.^2+F2.^2));
    lh=quiver3(XXplot,ZZplot,10*ones(size(XXplot)),F1F./modl,F2F./modl,zeros(size(XXplot)));
    
    set(lh,'linewidth',1.0);
    set(lh,'color',[0 0 0]);
    
    drawnow;
    disp(['time: ',num2str(time(iter+1)/86400/365.25)])
    pause;
end
