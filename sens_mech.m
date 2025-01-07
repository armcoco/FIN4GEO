clear all
close all
clc

FontSize=16;
% MarkerSize='default';
LineWidth='default';
MarkerFaceColor=[1 1 1];
set(0,'defaultaxesfontsize',FontSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

load sensitivity_analysis.mat
LW=1.01;
MS=14;
V_maxFB(10)=0.25963844;
figure
hold on
S1=0:0.2:3;
vvv=ones(size(S1));
vvv(1)=1e20;

plot(S1,log10(V_maxC(1:16)/V_maxC(1)),'d','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxC(1:16)/V_maxC(1)),1,vvv);
plot([0 3],[0 cc(1)*3],'-','lineWidth',LW)
cc(1)
plot(S1,log10(V_maxFA(1:16)/V_maxFA(1)),'o','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxFA(1:16)/V_maxFA(1)),1,vvv);
plot([0 3],[0 cc(1)*3],'--','lineWidth',LW)
cc(1)
plot(S1,log10(V_maxFB(1:16)/V_maxFB(1)),'*','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxFB(1:16)/V_maxFB(1)),1,vvv);
plot([0 3],[0 cc(1)*3],':','lineWidth',LW)
cc(1)
plot(S1,log10(V_maxC(16+(1:16))/V_maxC(17)),'dr','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxC(16+(1:16))/V_maxC(17)),1,vvv);
plot([0 3],[0 cc(1)*3],'-r','lineWidth',LW)
cc(1)
plot(S1,log10(V_maxFA(16+(1:16))/V_maxFA(17)),'or','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxFA(16+(1:16))/V_maxFA(17)),1,vvv);
plot([0 3],[0 cc(1)*3],'--r','lineWidth',LW)
cc(1)
plot(S1,log10(V_maxFB(16+(1:16))/V_maxFB(17)),'*r','lineWidth',LW,'markerSize',MS)
cc=wpolyfit(S1,log10(V_maxFB(16+(1:16))/V_maxFB(17)),1,vvv);
plot([0 3],[0 cc(1)*3],':r','lineWidth',LW)
cc(1)
legend('t=3 years - center of the model   ','linear bestfit at t=3 years - center of the model','t=3 years - Fault A','linear bestfit at t=3 years - Fault A', 't=3 years - Fault B','linear bestfit at t=3 years - Fault B',...
    't=100 years - center of the model','linear bestfit at t=100 years - center of the model','t=100 years - Fault A','linear bestfit at t=100 years - Fault A', 't=100 years - Fault B','linear bestfit at t=100 years - Fault B')
% plot([0,3],[0,0],'-k','lineWidth',LW)
box on
% xlabel('s_1 = \log_{10} \left( \bar{\mu} / \mu_c \right)')
% ylabel('s_2 = \log_{10} \left( v / v_0 \right)')
xlabel('s_1 = log_{10} ( \mu / \mu_c ) ','fontsize',20)
ylabel('s_2 = log_{10} ( v / v_0 )','fontsize',20)

return;

figure
hold on
plot(0:0.2:3,(V_maxC(1:16)/V_maxC(1)),'lineWidth',LW)
plot(0:0.2:3,(V_maxFA(1:16)/V_maxFA(1)),'--','lineWidth',LW)
plot(0:0.2:3,(V_maxFB(1:16)/V_maxFB(1)),':','lineWidth',LW)
plot(0:0.2:3,(V_maxC(16+(1:16))/V_maxC(17)),'r','lineWidth',LW)
plot(0:0.2:3,(V_maxFA(16+(1:16))/V_maxFA(17)),'--r','lineWidth',LW)
plot(0:0.2:3,(V_maxFB(16+(1:16))/V_maxFB(17)),':r','lineWidth',LW)
plot([0,3],[0,0],'-k','lineWidth',LW)
box on