clear all
close all

ColorSet=[1 0 0; 1 0.4 0; 0.8 0.8 0; 0 0.5 0; 0 0 1; 1 0 1];
LineStyles={'-','-','-','-'};
% LineStyles={'-','--','-.','.'};
Markers={'d','*','o','v'};

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];
Fault1P12=[Fault1P1(1)+(-1500-Fault1P1(2))/(Fault1P2(2)-Fault1P1(2))*(Fault1P2(1)-Fault1P1(1)) -1500];
Fault2P12=[Fault2P1(1)+(-1500-Fault2P1(2))/(Fault2P2(2)-Fault2P1(2))*(Fault2P2(1)-Fault2P1(1)) -1500];

% load Faults_paper_3Km_ScenarioA_x1.mat
% load Faults_paper_new_ScenarioA_x1.mat
% load Faults_paper_new_ScenarioA_x1_2yrs.mat
load Faults_paper_new_ScenarioA_x1_2yrs_4.mat
% load Faults_paper_new_ScenarioC_Gx0p1_x1.mat
% load Faults_paper_new_ScenarioC_Gx0p1_x1_516.mat
% load Faults_paper_new_ScenarioC_Gx10_x1_516.mat
% load Paoletta.mat
% load Faults_ScenarioC_x1_Trasatti.mat;
% Ntime=39;
% load Fault_inj_mass_constant_x2.mat
ttt=tt(selected_times+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
V_axis_onlyP(togli)=[];
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--r','lineWidth',2);
hold on
ttt=tt(selected_times+1);
ttt(togli)=[];
V_axis(togli)=[];
ttt=[0;ttt];
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,'--b','lineWidth',2)
plot(ttt',V_axis,'-k','lineWidth',2);
xlabel('time [years]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
legend('v_P (displacement due to pore pressure)','v_T (displacement due to thermal effects)','v=v_T+v_P (total displacement)')
% legend('v_P','v_T','v=v_T+v_P')
% title('UNREST \times 1','FontSize',20)

Vtime_onlyT=Vtime-Vtime_onlyP;
for iter=1:length(selected_times)
    perc(:,iter)=Vtime_onlyT(:,iter)./Vtime(:,iter);
    time_=round(tt(selected_times(iter)+1)*100)/100;
    if time_==1
        legend_str{iter}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter}=['time: ',num2str(time_),' years'];
    end
    
end

V_axis_onlyT=V_axis-V_axis_onlyP;
ToPlot={Vtime_onlyT,Vtime_onlyP,Vtime,Utime};
for j=1:length(ToPlot)
    ToPlot{j}=[2*ToPlot{j}(1,:)-1*ToPlot{j}(2,:);ToPlot{j}];
%     ToPlot{j}=interp1(Xusup,ToPlot{j},[0,Xusup]);
end
    
% ToPlot={[V_axis_onlyT(2:end);Vtime_onlyT],[V_axis_onlyP(2:end);Vtime_onlyP],[V_axis(2:end);Vtime],[U_axis;Utime]};
titles={'v_T (vertical displacement due to thermal effects) [m]',...
    'v_P (vertical displacement due to pore pressure) [m]',...
    'v=v_T+v_P (total vertical displacement) [m]',...
    'u (total horizontal displacement) [m]'};
ylabels={'vertical displacement [m]',...
    'vertical displacement [m]',...
    'vertical displacement [m]',...
    'horizontal displacement [m]'};
filesave={'defo_vT_ScenarioC_x1.eps',...
    'defo_vP_ScenarioC_x1.eps',...
    'defo_vPT_ScenarioC_x1.eps',...
    'defo_uPT_ScenarioC_x1.eps'};
plotABCD='DCAB';

selected_plots=[2,12,40,60]; % for Scenarios A, B and C
% selected_plots=[1,2,3,4]; % for Scenario Gx0p1
% ylims=[[-0.02 0.16];[-0.02 0.22];[-0.02 0.22];[-0.02 0.16]]; %for Scenarios A and B
ylims=[[-0.02 0.24];[-0.05 0.3];[-0.05 0.3];[-0.15 0.25]]; %for Scenario C
% ylims=[[-0.02 0.2];[-0.05 0.5];[-0.05 0.5];[-0.25 0.25]]; % for Scenario Gx0p1
% ylims=[[-0.02 0.32];[-0.05 0.25];[-0.05 0.4];[-0.25 0.35]]; % for Scenario Gx10

Rect=[0.1 0.1 0.8 0.8];
AxisPos = moPlotPos2(2,2, Rect,0.12,0.22);

figure1=figure;
plotPos=[3,1,2,4];
for plots=1:4
%     HH2=plot(Xusup(19:20:end),ToPlot{plots}(20:20:end,selected_plots),'lineWidth',2);
    axes('Position',AxisPos(plotPos(plots),:));
    HH1=plot([0,Xusup],ToPlot{plots}(:,selected_plots),'lineWidth',2);
    hold on
    xlim([0 10000])
    if plots<3
        xlabel('radial distance [m]','FontSize',20)
    else
        set(gca,'xticklabel',[])
    end
%     ylabel(ylabels{plots},'FontSize',20);
    title(titles{plots},'FontSize',20)
    for i=1:length(HH1)
        set(HH1(i),'color',ColorSet(i,:))
%         set(HH2(i),'color',ColorSet(i,:))
        set(HH1(i),'linestyle',LineStyles{i})
%         set(HH2(i),'linestyle','.')
%         set(HH2(i),'marker',Markers{i})
    end
    plot([Fault1P12(1),Fault1P12(1)],[-1 1],'-r','lineWidth',1.0)
    plot([Fault1P2(1),Fault1P2(1)],[-1 1],'-r','lineWidth',1.0)
    plot([Fault2P12(1),Fault2P12(1)],[-1 1],'-r','lineWidth',1.0)
    plot([Fault2P2(1),Fault2P2(1)],[-1 1],'-r','lineWidth',1.0)
    plot([200,200],[-1 1],'-r','lineWidth',1.0)
    if plots==4
        legend(legend_str{selected_plots})
    end
    ylim(ylims(plots,:))
%     print(filesave{i},'eps');
    text(3700,ylims(plots,2)-0.08*(ylims(plots,2)-ylims(plots,1)),['\color{black} \fontsize{48} ',plotABCD(plots)])
%     if plots>3
%        set(gca,'xticklabel',[]);
%     else
%         xlabel('radial distance [m]','FontSize',20);
%     end
%     saveas(figure1,filesave{plots},'epsc2');
end

return;

figure
HH=plot(Xusup,Vtime_onlyT(:,selected_plots),'lineWidth',2);
xlim([0 10000])
% legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('v_T (vertical displacement due to thermal effects)','FontSize',20)
hold on
for i=1:length(HH)
    set(HH(i),'color',ColorSet(i,:))
    set(HH(i),'linestyle',LineStyles{i})
    set(HH(i),'marker','d')
end
plot([Fault1P12(1),Fault1P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault1P2(1),Fault1P2(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P12(1),Fault2P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P2(1),Fault2P2(1)],[-1 1],'-r','lineWidth',1.0)
ylim([-0.02 0.16])

figure
HH1=plot(Xusup,Vtime_onlyP(:,selected_plots),'lineWidth',2);
hold on
HH2=plot(Xusup(10:10:end),Vtime_onlyP(10:10:end,selected_plots),'lineWidth',2);
xlim([0 10000])
% legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('v_P (vertical displacement due to pore pressure)','FontSize',20)
for i=1:length(HH)
    set(HH1(i),'color',ColorSet(i,:))
    set(HH2(i),'color',ColorSet(i,:))
    set(HH1(i),'linestyle',LineStyles{i})
    set(HH2(i),'linestyle','.')
    set(HH2(i),'marker',Markers{i})
end
plot([Fault1P12(1),Fault1P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault1P2(1),Fault1P2(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P12(1),Fault2P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P2(1),Fault2P2(1)],[-1 1],'-r','lineWidth',1.0)
ylim([-0.02 0.22])

figure
HH=plot(Xusup,Vtime(:,selected_plots),'lineWidth',2);
xlim([0 10000])
% legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('v=v_T+v_P (total vertical displacement)','FontSize',20)
hold on
for i=1:length(HH)
    set(HH(i),'color',ColorSet(i,:))
    set(HH(i),'linestyle',LineStyles{i})
    set(HH(i),'marker','d')
end
plot([Fault1P12(1),Fault1P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault1P2(1),Fault1P2(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P12(1),Fault2P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P2(1),Fault2P2(1)],[-1 1],'-r','lineWidth',1.0)
ylim([-0.02 0.22])

figure
HH=plot(Xusup,Utime(:,selected_plots),'lineWidth',2);
xlim([0 10000])
legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('horizontal displacement [m]','FontSize',20);
title('u (total horizontal displacement)','FontSize',20)
hold on
for i=1:length(HH)
    set(HH(i),'color',ColorSet(i,:))
    set(HH(i),'linestyle',LineStyles{i})
    set(HH(i),'marker','d')
end
plot([Fault1P12(1),Fault1P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault1P2(1),Fault1P2(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P12(1),Fault2P12(1)],[-1 1],'-r','lineWidth',1.0)
plot([Fault2P2(1),Fault2P2(1)],[-1 1],'-r','lineWidth',1.0)
ylim([-0.02 0.16])
