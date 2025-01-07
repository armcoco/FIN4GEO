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


% load Grav_3Km_ScenarioA_x1.mat
load Grav_paper_new_ScenarioA_x0p5.mat
dG_axis=dGtime(1,:);
ttt=tt(1:Ntime+1);
% togli=[8 10 13 14 16 17 20];
togli=[];
ttt(togli)=[];
dG_axis(togli)=[];
ttt=[0;ttt];
dG_axis=[0,dG_axis];
load Grav_paper_new_ScenarioA_x1.mat
dG_axis(2,:)=[0,dGtime(1,:)];
load Grav_paper_new_ScenarioA_x2.mat
dG_axis(3,:)=[0,dGtime(1,:)];
load Grav_paper_new_ScenarioA_x3.mat
dG_axis(4,:)=[0,dGtime(1,:)];
% HH2=plot(ttt([8,53:2:end])',-dG_axis(:,[8,53:2:end]),'-','lineWidth',2);
% HH2=plot(ttt([53:2:end])',-dG_axis(:,[53:2:end]),'-','lineWidth',2);
hold on
HH1=plot(ttt',-dG_axis,'-','lineWidth',2);
xlabel('time [years]','FontSize',20)
ylabel('gravity changes [\muGal]','FontSize',20);
legend('unrest \times 0.5','unrest \times 1','unrest \times 2','unrest \times 3');
for i=1:length(HH1)
    set(HH1(i),'color',ColorSet(i,:))
%     set(HH2(i),'color',ColorSet(i,:))
    set(HH1(i),'linestyle',LineStyles{i})
%     set(HH2(i),'linestyle','.')
%     set(HH2(i),'marker',Markers{i})
end
plot([0 2.5 2.5 0],[-350 -350 50 50],'k','linewidth',1.0)
% text(50,500-0.08*3500,['\color{black} \fontsize{48} B'])
text(0,-3000+0.1*3500,['\color{black} \fontsize{48} B'])
box on
% return;

% xlim([0 2.5]);

% HH3=plot(ttt(7)',-dG_axis(:,7),'-','lineWidth',2);
% for i=1:length(HH3)
%     set(HH3(i),'color',ColorSet(i,:))
%     set(HH3(i),'marker',Markers{i})
% end
% text(1.25,50-0.08*400,['\color{black} \fontsize{48} A'])
text(0,-350+0.1*400,['\color{black} \fontsize{48} A'])
legend off
box on
return;

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
