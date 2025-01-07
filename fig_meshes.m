clear all
close all

% LinearElasticity_header;

scale=1e3;

FontSize=16;
FS=20;
LW=3;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

load Xu.txt 
load Yu.txt

Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];

subplot(1,2,2)
plot(Xu/scale,Yu/scale,'b.')
hold on
xlim([0,50])
ylim([-50 0])

plot([Fault1P1(1),Fault1P2(1)]/scale,[Fault1P1(2),Fault1P2(2)]/scale,'-r','lineWidth',2)
plot([Fault2P1(1),Fault2P2(1)]/scale,[Fault2P1(2),Fault2P2(2)]/scale,'-r','lineWidth',2)

xlabel('radial distance [km]','fontSize',FS)
ylabel('z [km]','fontSize',FS)
title('mesh for the mechanical model ','fontSize',FS)

plot([0 10 10 0],[-1.5 -1.5 0 0],'color',[0.85 0.85 0],'lineWidth',2)
text(35,-35,['\color{black} \fontsize{48} B'])

%%%

subplot(1,2,1)
load xcenter.txt
load ycenter.txt
load zcenter.txt

plot(xcenter/scale,zcenter/scale,'b.')
hold on
plot([Fault1P1(1),Fault1P2(1)]/scale,[Fault1P1(2),Fault1P2(2)]/scale,'-r','lineWidth',2)
plot([Fault2P1(1),Fault2P2(1)]/scale,[Fault2P1(2),Fault2P2(2)]/scale,'-r','lineWidth',2)
xlim([0 10])
ylim([-1.5 0])
xlabel('radial distance [km]','fontSize',FS)
ylabel('z [km]','fontSize',FS)
title('mesh for the hydrological model ','fontSize',FS)
xlim([0,10])
ylim([-10 0])
plot([0 10 10 0],[-1.5 -1.5 0 0],'color',[0.85 0.85 0],'lineWidth',5)
text(7,-7,['\color{black} \fontsize{48} A'])
