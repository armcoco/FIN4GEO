close all

FA=0.5*(Fault1P12(1)+Fault1P2(1));
FB=0.5*(Fault2P12(1)+Fault2P2(1));

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

injcase={'0p5','1','2','3'};
for injs=1:4
    load(['Faults_paper_new_ScenarioB_x',num2str(injcase{injs}),'.mat'])
    ttt=tt(selected_times+1);
    ttt=[0;ttt];
    V_axis_onlyP=interp1(Xusup,Vtime_onlyP,FA);
    V_axis_onlyP=[0,V_axis_onlyP];
    plot(ttt',V_axis_onlyP,'-r','lineWidth',2);
    hold on
    V_axis=interp1(Xusup,Vtime,FA);
    V_axis=[0,V_axis];
    plot(ttt',V_axis-V_axis_onlyP,'-b','lineWidth',2)
    plot(ttt',V_axis,'-k','lineWidth',2);
    % ylim([-5e-4 0.1])
    xlabel('time [years]','FontSize',20)
    ylabel('vertical displacement [m]','FontSize',20);
    legend('v_P (displacement due to pore pressure)','v_T (displacement due to thermal effects)','v=v_T+v_P (total displacement)')
    % legend('v_P','v_T','v=v_T+v_P')
    % title('UNREST \times 1','FontSize',20)
end


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





