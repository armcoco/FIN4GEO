close all

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

load Faults_paper_3Km_ScenarioA_x1.mat
ttt=tt(selected_times+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
V_axis_onlyP(togli)=[];
ttt=[0;ttt];
V_axis_onlyP=[0,V_axis_onlyP];
ttt=tt(selected_times+1);
ttt(togli)=[];
V_axis(togli)=[];
ttt=[0;ttt];
V_axis=[0,V_axis];

Vtime_onlyT=Vtime-Vtime_onlyP;
for iter=1:length(selected_times)
    time_=round(tt(selected_times(iter)+1)*100)/100;
    if time_==1
        legend_str{iter}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter}=['time: ',num2str(time_),' years'];
    end
    
end

selected_plots=[1,3,5];

figure(1)
hold on
plot(Xusup,Vtime(:,selected_plots),'lineWidth',2)
xlim([0 10000])
% legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('vertical displacement [m]','FontSize',20);
title('v=v_T+v_P (total vertical displacement)','FontSize',20)

figure(2)
hold on
plot(Xusup,Utime(:,selected_plots),'lineWidth',2)
xlim([0 10000])
legend(legend_str{selected_plots})
xlabel('radial distance [m]','FontSize',20)
ylabel('horizontal displacement [m]','FontSize',20);
title('u (total horizontal displacement)','FontSize',20)

load Faults_paper_3Km_ScenarioA_x1_noFaults.mat

figure(1)
plot(Xusup,Vtime(:,selected_plots),'--','lineWidth',2)
xlim([0 10000])
figure(2)
plot(Xusup,Utime(:,selected_plots),'--','lineWidth',2)
xlim([0 10000])
