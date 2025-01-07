clear all
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

% pathSol='Faults_paper_new/ScenarioC_2yrs/unrest_x1/';
% pathSol='Faults_paper_new/ScenarioC/unrest_x1/';
% pathSol='Faults_paper_3Km/ScenarioA/unrest_x1/';
% pathSol='three_layers/ScenarioA/unrest_x1/';
% pathSol='EOS2_only_water/unrest_x1_6100/';
pathSol='EOS2_only_water/unrest_x1/';


Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];


Sol=load([pathSol,'Sol_',num2str(steady)]);
X=Sol(:,1);
Z=Sol(:,3);
T0=Sol(:,4);
P0=Sol(:,5);
a=find(Z==Z(1));
l=length(a);
m=length(Z)/l;
TT0=reshape(T0,l,m);
PP0=reshape(P0,l,m);
Xcenter=load([pathSol,'xcenter.txt']);
Ycenter=load([pathSol,'ycenter.txt']);
Zcenter=load([pathSol,'zcenter.txt']);
[porosity,rho_rock]=por_and_rho(Xcenter,Ycenter,Zcenter,Fault1P1,Fault1P2,Fault2P1,Fault2P2);
Sg0=Sol(:,6);
Dw0=Sol(:,7);
Dg0=Sol(:,8);
Rho0=porosity.*(Dw0.*(1-Sg0)+Dg0.*Sg0);
RRho0=reshape(Rho0,l,m);
SSGG0=reshape(Sg0,l,m);

% time=load('faults_3000_rainfall_only/unrest_and_quiet_rainfall_constant/time');
time=load([pathSol,'time']);
for iter=0:length(time)-1 %iter=[1,11,39,59] %
    time_=time(iter+1)/86400/365.25;
    clf
    Sol=load([pathSol,'Sol_',num2str(iter)]);
    Conn=load([pathSol,'Conn_',num2str(iter)]);
    X=Sol(:,1);
    Z=Sol(:,3);
    T=Sol(:,4);
    
    P=Sol(:,5);
    SG=Sol(:,6);
    Sg=Sol(:,6);
    Dw=Sol(:,7);
    Dg=Sol(:,8);
    Rho=porosity.*(Dw.*(1-Sg)+Dg.*Sg);

    a=find(Z==Z(1));
    l=length(a);
    m=length(Z)/l;
    XX=reshape(X,l,m);
    ZZ=reshape(Z,l,m);
    TT=reshape(T,l,m);
    PP=reshape(P,l,m);
    SSGG=reshape(SG,l,m);
    RRho=reshape(Rho,l,m);
    
    subplot(3,1,1)
    [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
%     clabel(c,h);
%     caxis([-0.5 30])
    caxis([0 2.5])
    %     xlim([0 1000])
    colorbar
    title(['pore pressure changes P-P_0 [MPa],   time: ',num2str(round(time_*100)/100),' years'],'FontSize',16)
    
    subplot(3,1,2)
    [c,h]=contourf(XX,ZZ,TT-TT0);
%     clabel(c,h);
%     caxis([0 360])
    caxis([-15 100])
    caxis([-15 170])
    caxis([0 5])
    %     xlim([0 1000])
    colorbar
    title('temperature changes T-T_0 [C]','FontSize',16)
    hold on
%     plottaGradT;
    
    subplot(3,1,3)
%     [c,h]=contourf(XX,ZZ,-(RRho-RRho0));
    [c,h]=contourf(XX,ZZ,SSGG-SSGG0);
%     clabel(c,h);
%     caxis([0 360])
    caxis([-15 200])
    caxis([-15 320])
    caxis([0 1])
    %     xlim([0 1000])
    colorbar
%     title('density changes -(\rho-\rho_0) [Kg/m^3]','FontSize',16)
    title('gas saturation changes S_g-S_{g0}','FontSize',16)
    
    drawnow;
    disp(['time: ',num2str((time(iter+1)/86400/365.25)),' years'])
    pause;
end
