clear all
close all

FontSize=14;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


% load Grav_3Km_ScenarioA_x1.mat
% load Grav_paper_new_ScenarioC_x1.mat
% load Faults_paper_new_ScenarioC_x1.mat
load Faults_paper_new_ScenarioA_x1_516.mat
load Grav_paper_new_ScenarioA_x1_516.mat
% load Paola.mat
ttt=tt(1:Ntime+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
ttt=[0;ttt];

for iter=0:Ntime
    time_=round(tt(iter+1)*100)/100;
    if time_==1
        legend_str{iter+1}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter+1}=['time: ',num2str(time_),' years'];
    end    
end

selected_plots=[1,2,4,5];
selected_plots=1:6; %4:6;
hFig=figure;
co = get(gca,'ColorOrder');

plot(Xusup(2:end-1),-dGtime(2:end-1,selected_plots)./Vtime(:,selected_plots),'lineWidth',2);
% [AX,H1,H2]=plotyy(Xusup(2:end-1),-dGtime(2:end-1,selected_plots),Xusup(2:end-1),Vtime(:,selected_plots),'plot','plot');
set(gcf, 'Position', [100 1000 400 400]);
set(gca, 'Position', [0.2 0.12 0.6 0.6]);
xlim([0,8000]);
% ylim([-20000 20000]);
% set(get(AX(1),'Ylabel'),'String','gravity changes [\muGal] (cont. line)','FontSize',FontSize*1.2);
xlabel('radial distance [m]','FontSize',FontSize*1.2)
ylabel('gravity gradient [\muGal/m]','FontSize',FontSize*1.2);
legend(legend_str{selected_times(selected_plots)+1})
