close all

ColorSet=[1 0 0; 1 0.4 0; 0.8 0.8 0; 0 0.5 0; 0 0 1; 1 0 1];

Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];
Fault1P12=[Fault1P1(1)+(-1500-Fault1P1(2))/(Fault1P2(2)-Fault1P1(2))*(Fault1P2(1)-Fault1P1(1)) -1500];
Fault2P12=[Fault2P1(1)+(-1500-Fault2P1(2))/(Fault2P2(2)-Fault2P1(2))*(Fault2P2(1)-Fault2P1(1)) -1500];

FA=0.5*(Fault1P12(1)+Fault1P2(1));
FB=0.5*(Fault2P12(1)+Fault2P2(1));

FontSize=16;
LW=3;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


load Faults_paper_new_ScenarioA_x0p5.mat
ttt=tt(selected_times+1);
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--','Color',ColorSet(1,:),'lineWidth',LW);
hold on
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,':','Color',ColorSet(1,:),'lineWidth',LW)
H1=plot(ttt',V_axis,'-','Color',ColorSet(1,:),'lineWidth',LW);

load Faults_paper_new_ScenarioA_x1.mat
ttt=tt(selected_times+1);
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--','Color',ColorSet(2,:),'lineWidth',LW);
hold on
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,':','Color',ColorSet(2,:),'lineWidth',LW)
H2=plot(ttt',V_axis,'-','Color',ColorSet(2,:),'lineWidth',LW);
% ylim([-5e-4 0.1])
xlabel('time [years]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
% legend('v_P (displacement due to pore pressure)','v_T (displacement due to thermal effects)','v=v_T+v_P (total displacement)')
% legend('v_P','v_T','v=v_T+v_P')
% title('UNREST \times 1','FontSize',20)

load Faults_paper_new_ScenarioA_x2.mat
ttt=tt(selected_times+1);
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--','Color',ColorSet(3,:),'lineWidth',LW);
hold on
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,':','Color',ColorSet(3,:),'lineWidth',LW)
H3=plot(ttt',V_axis,'-','Color',ColorSet(3,:),'lineWidth',LW);

load Faults_paper_new_ScenarioA_x3.mat
ttt=tt(selected_times+1);
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--','Color',ColorSet(4,:),'lineWidth',LW);
hold on
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,':','Color',ColorSet(4,:),'lineWidth',LW)
H4=plot(ttt',V_axis,'-','Color',ColorSet(4,:),'lineWidth',LW);

legend_str={'unrest \times 0.5','unrest \times 1','unrest \times 2','unrest \times 3'};
legend([H1,H2,H3,H4],legend_str)

return; 

text(80,0.647,'\color{black} \fontsize{18} Unrest x3')
text(80,0.455,'\color{black} \fontsize{18} Unrest x2')
text(80,0.37,'\color{blue} \fontsize{18} Unrest x3')
text(80,0.295,'\color{red} \fontsize{18} Unrest x3')
text(80,0.25,'\color{blue} \fontsize{18} Unrest x2')
text(80,0.218,'\color{black} \fontsize{18} Unrest x1')
text(80,0.165,'\color{red} \fontsize{18} Unrest x2')
text(80,0.14,'\color{blue} \fontsize{18} Unrest x1')
text(80,0.09,'\color{red} \fontsize{18} Unrest x1')
text(80,0.06,'\color{black} \fontsize{18} Unrest x0.5')
text(80,0.04,'\color{red} \fontsize{18} Unrest x0.5')
text(80,0.01,'\color{blue} \fontsize{18} Unrest x0.5')





