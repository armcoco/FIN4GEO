clear all
close all

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


% load Faults_3Km_paper_ScenarioA_x1.mat
% load Grav_3Km_paper_ScenarioA_x1.mat
load Faults_paper_new_ScenarioB_x1_516.mat
load Grav_paper_new_ScenarioB_x1_516.mat
% load Paola_defoC2.mat
% load Paola_gravC2.mat
dG_axis=dGtime(1,:);
ttt=tt(1:Ntime+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
dG_axis(togli)=[];
ttt=[0;ttt];
dG_axis=[0,dG_axis];
% % % plot(ttt',dG_axis,'-b','lineWidth',2);
% % % hold on
% % % xlabel('time [years]','FontSize',20)
% % % ylabel('gravity changes [\muGal]','FontSize',20);


for iter=0:Ntime
    time_=round(tt(iter+1)*100)/100;
    if time_==1
        legend_str{iter+1}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter+1}=['time: ',num2str(time_),' years'];
    end    
end

% selected_times=[2,12,20,40,52];
% selected_times=1:5;
% selected_times=[40];
selected_plots=1:6;
figure
plot(Xusup(2:end-1),dGtime(2:end-1,selected_plots)./Vtime(:,selected_plots),'lineWidth',2)
xlim([0 10000])
legend(legend_str{selected_times+1})
xlabel('radial distance [m]','FontSize',20)
ylabel('gravity gradient [\muGal/m]','FontSize',20);
