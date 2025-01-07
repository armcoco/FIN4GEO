clear all
close all

ColorSet=[1 0 0; 1 0.4 0; 0.8 0.8 0; 0 0.5 0; 0 0 1; 1 0 1];

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
selected_plots=4:6;
hFig=figure;
co = get(gca,'ColorOrder');
[AX,H1,H2]=plotyy(Xusup(2:end-1),-dGtime(2:end-1,selected_plots),Xusup(2:end-1),Vtime(:,selected_plots),'plot','plot');
set(gcf, 'Position', [100 1000 400 400]);
set(gca, 'Position', [0.2 0.12 0.6 0.6]);
set(AX,'xlim',[0,8000]);
set(AX(1),'ylim',[-300 300],'ytick',[-300:100:300]);
% set(AX(2),'ylim',[-0.35 0.35],'ytick',[-0.35:0.07:0.35]);
set(AX(2),'ylim',[-0.3 0.3],'ytick',[-0.3:0.1:0.3]);
% set(get(AX(1),'Ylabel'),'String','gravity changes [\muGal] (continuous line)','FontSize',FontSize);
% set(get(AX(2),'Ylabel'),'String','vertical displacement [m] (dashed line)','FontSize',FontSize)
set(get(AX(1),'Ylabel'),'String','gravity changes [\muGal] (cont. line)','FontSize',FontSize*1.2);
set(get(AX(2),'Ylabel'),'String','vertical displ. [m] (dashed line)','FontSize',FontSize*1.2)
xlabel('radial distance [m]','FontSize',FontSize*1.2)
% title('qqq')
set(H1,'LineStyle','-')
set(H2,'LineStyle','--')
set(H1,'lineWidth',2)
set(H2,'lineWidth',2)
for i=1:length(selected_plots)
    set(H1(i),'color',co(i+3,:))
    set(H2(i),'color',co(i+3,:))
end
legend(legend_str{selected_times(selected_plots)+1})
