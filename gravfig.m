clear all
close all

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


% load Grav_3Km_ScenarioA_x1.mat
load Grav_paper_new_ScenarioA_x1.mat
% load Paola.mat
dG_axis=dGtime(1,:);
ttt=tt(1:Ntime+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
dG_axis(togli)=[];
ttt=[0;ttt];
dG_axis=[0,dG_axis];
plot(ttt',dG_axis,'-b','lineWidth',2);
hold on
xlabel('time [years]','FontSize',20)
ylabel('gravity changes [\muGal]','FontSize',20);


for iter=0:Ntime
    time_=round(tt(iter+1)*100)/100;
    if time_==1
        legend_str{iter+1}=['time: ',num2str(time_),' year'];
    else
        legend_str{iter+1}=['time: ',num2str(time_),' years'];
    end    
end

selected_times=[2,12,40,60];
figure
plot(Xusup(2:end-1),dGtime(2:end-1,selected_times),'lineWidth',2)
xlim([0 10000])
legend(legend_str{selected_times})
xlabel('radial distance [m]','FontSize',20)
ylabel('gravity changes [\muGal]','FontSize',20);
