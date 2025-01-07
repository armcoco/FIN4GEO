clear all
close all

ColorSet=[1 0 0; 1 0.4 0; 0.8 0.8 0; 0 0.5 0; 0 0 1; 0.7 0 0.7];

FontSize=12;
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

Scenario_letter='ABC';
Scenario_sym={'I','II','III'};
plotABCDEFGHI='GDAHEBIFC';

load Faults_paper_new_ScenarioA_x1_516.mat
load Grav_paper_new_ScenarioA_x1_516.mat
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

nRows=3;
nCols=3;
Rect=[0.08,0.06,0.84,0.88];
AxisPos = moPlotPos2(nCols,nRows,Rect);
figure
iterfig=1;
for scenario=1:3
    load(['Faults_paper_new_Scenario',Scenario_letter(scenario),'_x1_516.mat'])
    load(['Grav_paper_new_Scenario',Scenario_letter(scenario),'_x1_516.mat'])
    for rowplot=1:nRows
        axes('Position',AxisPos(rowplot+(scenario-1)*nRows,:));
        if rowplot<3
            selected_plots=(1:3)+3*(rowplot==1);
            [AX,H1,H2]=plotyy(Xusup(2:end-1),-dGtime(2:end-1,selected_plots),Xusup(2:end-1),Vtime(:,selected_plots),'plot','plot');
        else
            selected_plots=1:6;
            if scenario==1
                dGtime(Xusup>3500,selected_plots)=NaN;
                Vtime(Xusup(2:end-1)>3500,selected_plots)=NaN;
                Xusup(Xusup>3500)=NaN;
            end
            H1=plot(Xusup(2:end-1),-dGtime(2:end-1,selected_plots)./Vtime(:,selected_plots),'lineWidth',2);
            H2=H1;
            AX=gca;
        end
        hold on
        plot([Fault1P12(1),Fault1P12(1)],[-1 1]*1e5,'-r','lineWidth',1.0)
        plot([Fault1P2(1),Fault1P2(1)],[-1 1]*1e5,'-r','lineWidth',1.0)
        plot([Fault2P12(1),Fault2P12(1)],[-1 1]*1e5,'-r','lineWidth',1.0)
        plot([Fault2P2(1),Fault2P2(1)],[-1 1]*1e5,'-r','lineWidth',1.0)
        plot([200,200],[-1 1]*1e5,'-r','lineWidth',1.0)
        set(AX,'xlim',[0,8000],'xtick',[0:1000:8000]);
        if rowplot==1
            set(AX(1),'ylim',[-1600 1600],'ytick',[-1600:400:1600]);
            set(AX(1),'yticklabel',{'','-1200','','-400','','400','','1200',''})
            set(AX(2),'ylim',[-0.4 0.4],'ytick',[-0.4:0.1:0.4]);
            set(AX(2),'yticklabel',{'','-0.3','','-0.1','','0.1','','0.3',''})
            xlabel('radial distance [km]','FontSize',FontSize);
            set(AX,'xticklabel',{'','1','','3','','5','','7',''})
        elseif rowplot==2
            set(AX(1),'ylim',[-300 300],'ytick',[-300:100:300]);
            set(AX(1),'yticklabel',{'','  -200','','0','','  200',''})
            set(AX(2),'ylim',[-0.3 0.3],'ytick',[-0.3:0.1:0.3]);
            set(AX(2),'yticklabel',{'','-0.2','','0','','0.2',''})
            set(AX,'xticklabel',[])
        elseif rowplot==3
            set(AX,'ylim',[-20000 1000],'ytick',-18000:3000:0);
%             set(AX,'yticklabel',{'-1800','','','0','','  200',''})
            set(AX,'xticklabel',[])
        end
        if scenario==1
            if rowplot<3
                set(get(AX(1),'Ylabel'),'String','gravity changes [\muGal] (solid line)','FontSize',FontSize);
                set(AX(2),'yticklabel',[])
            else
                set(get(AX,'Ylabel'),'String','gravity gradient [\muGal/m]','FontSize',FontSize);
            end
        elseif scenario==2
            set(AX,'yticklabel',[])
        elseif scenario==3
            if rowplot<3
                set(AX(1),'yticklabel',[])
                set(get(AX(2),'Ylabel'),'String','vertical displ. [m] (dashed line)','FontSize',FontSize)
            else
                set(AX,'yticklabel',[])
            end
        end
        % title('qqq')
        set(H2,'LineStyle','--')
        set(H1,'LineStyle','-')
        set(H1,'lineWidth',2)
        set(H2,'lineWidth',2)
        for i=1:length(selected_plots)
            set(H1(i),'color',ColorSet(selected_plots(i),:))
            set(H2(i),'color',ColorSet(selected_plots(i),:))
        end
        if scenario==1
            legend(legend_str{selected_times(selected_plots)+1})
        end
        if rowplot==3
            title(['SCENARIO ',Scenario_sym{scenario}],'FontSize',18);
        end
        xl=xlim;
        yl=ylim;
        text(mean(xl),yl(1)+0.14*(yl(2)-yl(1)),['\color{black} \fontsize{36} ',plotABCDEFGHI(iterfig)]);
        iterfig=iterfig+1;
    end
end
