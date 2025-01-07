close all

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


load Faults_paper_ScenarioB_x1.mat
% Ntime=39;
% load Fault_inj_mass_constant_x2.mat 
ttt=tt(1:Ntime+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
V_axis_onlyP(togli)=[];
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
plot(ttt',V_axis_onlyP,'--r','lineWidth',2);
hold on
ttt=tt(1:Ntime+1);
ttt(togli)=[];
V_axis(togli)=[];
ttt=[0;ttt];
V_axis=[0,V_axis];
plot(ttt',V_axis-V_axis_onlyP,'-.b','lineWidth',2)
plot(ttt',V_axis,'-k','lineWidth',2);
xlabel('time [years]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
legend('u_P (displacement due to pore pressure)','u_T (displacement due to thermal effects)','u=u_T+u_P (total displacement)')


% return;


Vtime_onlyT=Vtime-Vtime_onlyP;
for iter=0:Ntime
    perc(:,iter+1)=Vtime_onlyT(:,iter+1)./Vtime(:,iter+1);
    time_=round(tt(iter+1)*100)/100;
    if time_==1
        legend_str{iter+1}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter+1}=['time: ',num2str(time_),' years'];
    end
    
end
% selected_times=[2,4,12,20,40,55,60];
selected_times=[2,4,12,20,40,60];
selected_times=[2,12,40,60];

figure
plot(Xusup,Vtime_onlyT(:,selected_times),'lineWidth',2)
xlim([0 10000])
% legend(legend_str{selected_times})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('u_T (vertical displacement due to thermal effects)','FontSize',20)

figure
plot(Xusup,Vtime_onlyP(:,selected_times),'lineWidth',2)
xlim([0 10000])
% legend(legend_str{selected_times})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('u_P (vertical displacement due to pore pressure)','FontSize',20)

figure
plot(Xusup,Utime(:,selected_times),'lineWidth',2)
xlim([0 10000])
legend(legend_str{selected_times})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('u=u_T+u_P (total vertical displacement)','FontSize',20)
