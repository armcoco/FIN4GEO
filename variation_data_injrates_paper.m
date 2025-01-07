clear all
close all

FontSize = 12;
FontSize2 = 16;
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

nScenarios=3;
% pathSol='Faults_paper_new/ScenarioC_2yrs/unrest_x1/';
pathSol=['Faults_paper_new/ScenarioA/unrest_x1/'];
pathSol3{1}='Faults_paper_new/ScenarioC/unrest_x1/';
pathSol3{2}='Faults_paper_new/ScenarioB/unrest_x1/';
pathSol3{3}='Faults_paper_new/ScenarioA/unrest_x1/';
pathSol3{1}='Faults_paper_new/ScenarioA/unrest_x0p5/';
pathSol3{2}='Faults_paper_new/ScenarioA/unrest_x2/';
pathSol3{3}='Faults_paper_new/ScenarioA/unrest_x3/';
% Scenario_letter={'III','II','I'};
Scenario_letter={'0.5','2','3'};
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

Rwidth=0.3;
Rheight=0.22;
Rwidth_off=0.01;
Rheight_off=0.02;
Rx=0.040+(0:2)*(Rwidth+Rwidth_off);
Ry=0.055+(0:2)*(Rheight+Rheight_off);
for jj=1:3
    for ii=1:3
        Rect{ii,jj} = [Rx(jj),Ry(ii),Rwidth,Rheight];
    end
end

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

for variable=1:3
    for scenario=1:nScenarios
        pathSol=pathSol3{scenario};
        nplot=1;
        AxisPos = moPlotPos(1, 4, Rect{scenario,variable});
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
            
            axes('Position', AxisPos(nplot, :));
            %     clabel(c,h);
            if variable ==1
                [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
                hold on
                caxis([0 15]) %pressure
            elseif variable ==2
                [c,h]=contourf(XX,ZZ,(TT-TT0));
                hold on
                caxis([-15 190]) %temperature
            elseif variable ==3
                [c,h]=contourf(XX,ZZ,(SSGG-SSGG0));
                hold on
                caxis([0 1]) %gas saturation
            end
            %     xlim([0 1000])
            %     colorbar
            disp(['time: ',num2str((time(iter+1)/86400/365.25)),' years'])
            plot([Fault1P1(1),Fault1P2(1)],[Fault1P1(2),Fault1P2(2)],'-w','lineWidth',1.2)
            plot([Fault2P1(1),Fault2P2(1)],[Fault2P1(2),Fault2P2(2)],'-w','lineWidth',1.2)
            % ax = gca;
            % set(gca, 'XTickLabelMode', 'Manual')
            %      set(gca, 'XTick', [1000:1000:5000])
            % ax.YTick = [-1500 -1000 -500 0];
            % ax.XTick = [1000];
            
            if nplot==1
%                 text(4000,-500,['\color{white} \fontsize{14} Scenario ',Scenario_letter{scenario}]);
                text(4000,-500,['\color{white} \fontsize{14} Unrest  \times ',Scenario_letter{scenario}]);
                if scenario==3
                    %             title(['pore pressure changes P-P_0 [MPa],   time: ',num2str(round(time_*100)/100),' years'],'FontSize',16)
                    if variable ==1
                        title(['Pore pressure changes   \Delta P = P-P_0 [MPa]'],'FontSize',FontSize2)
                    elseif variable ==2
                        title(['Temperature changes   \Delta T = T-T_0 [C]'],'FontSize',FontSize2)
                    elseif variable ==3
                        title(['Gas saturation changes   \Delta S_g = S_g-S_{g0}'],'FontSize',FontSize2)
                    end
                    Rect_temp=Rect{nScenarios,variable};
                    colorbar('location','northoutside','position',[Rect_temp(1), Rect_temp(2)+Rect_temp(4)+0.04, Rect_temp(3), 0.01]);
                    colormap default
                end
            end
            if nplot<4
                set(gca,'xticklabel',[])
            end
            if nplot==4 && scenario==1
                xlabel('radial distance [km]','FontSize',FontSize2)
                set(gca,'xticklabel',{'1','2','3','4','5','6','7','8','9'})
            else
                set(gca,'xticklabel',[])
            end
            if variable==1
                set(gca,'yticklabel',{'-1','-0.5','0'})
                if nplot==3
                    ylabel('       z [km]','FontSize',FontSize2)
                end
            else
                set(gca,'yticklabel',[])
            end
            text(7500,-500,['\color{white} \fontsize{14} t=',num2str(round(time_*100)/100),' years']);
            nplot=nplot+1;
        end
    end
    Rect_temp=Rect{nScenarios,variable};
    axes('Position', [Rect_temp(1), Rect_temp(2)+Rect_temp(4)+0.09, Rect_temp(3), Rwidth/4*0.9]);
    %     clabel(c,h);
    if variable ==1
        [c,h]=contourf(XX,ZZ,(PP0)/1e6);
        hold on
%         caxis([0 4]) %pressure
        title(['Initial condition: pore pressure P_0 [MPa]'],'FontSize',FontSize2)
    elseif variable ==2
        [c,h]=contourf(XX,ZZ,(TT0));
        hold on
%         caxis([-15 170]) %temperature
        title(['Initial condition: temperature T_0 [C]'],'FontSize',FontSize2)
        set(gca,'yticklabel',[])
    elseif variable ==3
        [c,h]=contourf(XX,ZZ,(SSGG0));
        hold on
        caxis([0 1]) %gas saturation
        title(['Initial condition: gas saturation S_{g0}'],'FontSize',FontSize2)
        set(gca,'yticklabel',[])
    end
                plot([Fault1P1(1),Fault1P2(1)],[Fault1P1(2),Fault1P2(2)],'-w','lineWidth',1.2)
            plot([Fault2P1(1),Fault2P2(1)],[Fault2P1(2),Fault2P2(2)],'-w','lineWidth',1.2)
    set(gca,'xticklabel',[])
    colorbar('location','northoutside','position',[Rect_temp(1), Rect_temp(2)+Rect_temp(4)+0.2, Rect_temp(3), 0.01]);
    colormap default
end
