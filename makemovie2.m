% clear all
close all

FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


load Rinaldi3.mat
steady=1000;
% unrest=10;
% endsim=48;

Sol=load(['Faults/unrest/Sol_',num2str(steady)]);
X=Sol(:,1);
Z=Sol(:,3);
T0=Sol(:,4);
P0=Sol(:,5);
a=find(Z==Z(1));
l=length(a);
m=length(Z)/l;
TT0=reshape(T0,l,m);
PP0=reshape(P0,l,m);


% time=load('faults_3000_rainfall_only/unrest_and_quiet_rainfall_constant/time');
time=load('Faults/unrest/time');
for iter=0:length(time)-1
    time_=time(iter+1)/86400/365.25;
    clf
    Sol=load(['Faults/unrest/Sol_',num2str(iter)]);
    Conn=load(['Faults/unrest/Conn_',num2str(iter)]);
    X=Sol(:,1);
    Z=Sol(:,3);
    T=Sol(:,4);
    
    P=Sol(:,5);
    SG=Sol(:,6);
    
    a=find(Z==Z(1));
    l=length(a);
    m=length(Z)/l;
    XX=reshape(X,l,m);
    ZZ=reshape(Z,l,m);
    TT=reshape(T,l,m);
    PP=reshape(P,l,m);
    SSGG=reshape(SG,l,m);
    
    subplot(3,1,1)
    [c,h]=contourf(XX,ZZ,PP0/1e6);
    %     clabel(c,h);
    %     caxis([-0.5 30])
    caxis([0 15])
        xlim([0 3000])
    colorbar
    title(['pressure at time=0 [MPa]'],'FontSize',16)
    
    subplot(3,1,2)
    [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
    %     clabel(c,h);
    %     caxis([0 360])
    caxis([-1 4])
        xlim([0 3000])
    colorbar
    title('pressure change: P-P_0 [MPa]','FontSize',16)
    hold on
    %     plottaGradT;
    
    subplot(3,1,3)
    ttt=tt(1:Ntime+1);
    togli=[8 10 13 14 16 17 20];
    ttt(togli)=[];
    V_axis_onlyP_new=V_axis_onlyP;
    V_axis_onlyP_new(togli)=[];
    ttt=[0;ttt];
    V_axis_onlyP_new=[0,V_axis_onlyP_new];
%     plot(ttt',V_axis_onlyP,'--r','lineWidth',2);
    hold on
    ttt=tt(1:Ntime+1);
    ttt(togli)=[];
    V_axis_new=V_axis;
    V_axis_new(togli)=[];
    ttt=[0;ttt];
    V_axis_new=[0,V_axis_new];
    plot(ttt',V_axis_onlyP_new,'-b','lineWidth',2)
    hold on
    plot(tt(iter+1)',V_axis_onlyP(iter+1),'or','lineWidth',2)
    ylim([-0.02 0.20])
    
    drawnow;
    disp(['time: ',num2str(time(iter+1)/86400/365.25)])
%     pause;
    movie(iter+1)=getframe(gcf);
end
