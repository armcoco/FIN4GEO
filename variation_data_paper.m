clear all
close all

FontSize = 20;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];

steady=1000;
% unrest=10;
% endsim=48;

% pathSol='Faults_paper_new/ScenarioC_2yrs/unrest_x1/';
pathSol='Faults_paper_new/ScenarioA/unrest_x1/';
% pathSol='Faults_paper_3Km/ScenarioA/unrest_x1/';
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
XX=reshape(X,l,m);
ZZ=reshape(Z,l,m);
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

% figure;
% Rect = [0.19, 0.07, 0.775, 0.845];
Rect = [0.12, 0.11, 0.775, 0.695];
AxisPos = moPlotPos(1, 4, Rect);
% for i = 1:16
%   axes('Position', AxisPos(i, :));
% end

% subplot(5,1,1)
% [c,h]=contourf(XX,ZZ,SSGG0/1e6);
% clabel(c,h);
% caxis([-0.3 2])
% xlim([0 1000])
% colorbar
% title('Pore pressure - initial condition [MPa]','FontSize',16)
% title('Temperature - initial condition [C]','FontSize',16)
% title('Gas saturation - initial condition','FontSize',16)
% ylabel('z [m]','FontSize',16)
% set(gca,'xticklabel',[])

% time=load('faults_3000_rainfall_only/unrest_and_quiet_rainfall_constant/time');
time=load([pathSol,'time']);
figure
nplot=1;
for iter=[1,11,39,59] %0:length(time)-1
    time_=time(iter+1)/86400/365.25;
%     clf
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
    
%     subplot(5,1,nplot+1)
%     PP(XX>4000)=PP0(XX>4000);
% [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
%     [c,h]=contourf(XX,ZZ,(TT-TT0));
axes('Position', AxisPos(nplot, :));
    [c,h]=contourf(XX,ZZ,(SSGG-SSGG0));
    hold on
%     clabel(c,h);
%     caxis([0 4]) %pressure
%     caxis([-15 170]) %temperature
    caxis([0 1]) %gas saturation
    %     xlim([0 1000])
%     colorbar
    disp(['time: ',num2str((time(iter+1)/86400/365.25)),' years'])
    ylabel('z [m]','FontSize',26)
    plot([Fault1P1(1),Fault1P2(1)],[Fault1P1(2),Fault1P2(2)],'--w','lineWidth',1.1)
plot([Fault2P1(1),Fault2P2(1)],[Fault2P1(2),Fault2P2(2)],'--w','lineWidth',1.1)
% ax = gca;
% set(gca, 'XTickLabelMode', 'Manual')
%      set(gca, 'XTick', [1000:1000:5000])
% ax.YTick = [-1500 -1000 -500 0];
% ax.XTick = [1000];

if nplot==1
%             title(['pore pressure changes P-P_0 [MPa],   time: ',num2str(round(time_*100)/100),' years'],'FontSize',16)
%             title(['Pore pressure changes P-P_0 [MPa]'],'FontSize',16)
%             title(['Temperature changes T-T_0 [C]'],'FontSize',16)
            title(['Gas saturation changes S_g-S_{g0}'],'FontSize',26)
            
%         colorbar('northoutside');
    end
    if nplot<4
        set(gca,'xticklabel',[])
    end
    if nplot==4
        xlabel('radial distance [m]','FontSize',26)
    end
    nplot=nplot+1;
end
