clear all
close all
clc

load Serapeo.txt
% load Grav_paper_new_ScenarioA_x0p5.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x0p5.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17_hybrid.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17_linear_decay.mat
load Grav_paper_new_extended_release_200m_ScenarioA_x1p15_molar_ratio_0p17_Inj2000m_hybrid.mat
dG_axis=dGtime(1,:);
dG_axis=[0,dG_axis];
% dG_axis=max(0,dG_axis);
% delay=1*1.75;
delay=1*2.00;
ttt=[0;tt];
Ser=interp1(Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,ttt);

plot(ttt,-dG_axis)
hold on
xlim([0 20])
% plot(ttt,Ser)
plot(Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,'-*')
df=Ser+dG_axis';
% plot(ttt,df,'k')
defo_magma=660*(1-exp(-0.23*ttt));
plot(ttt,defo_magma','k')
plot(ttt,defo_magma-dG_axis(:),'r')

legend('gravity changes due to hydrothermal effects','gravity residual - real data','gravity changes due to magma','gravity changes for model A')
title('Serapeo')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load LaPietra.txt
figure
dG_axis=dGtime(18,:);
dG_axis=[0,dG_axis];
% dG_axis=max(0,dG_axis);
delay=1*1.75;
delay=1*2.00;
ttt=[0;tt];
Solf=interp1(LaPietra(:,1)/365.25-delay,LaPietra(:,2)*1e3,ttt);

plot(ttt,-dG_axis)
hold on
xlim([0 20])
% plot(ttt,Ser)
plot(LaPietra(:,1)/365.25-delay,LaPietra(:,2)*1e3,'-*')
df=Solf+dG_axis';
% plot(ttt,df,'y')
defo_magma=130*(1-exp(-0.5*ttt));
plot(ttt,defo_magma','k')
plot(ttt,defo_magma-dG_axis(:),'r')

legend('gravity changes due to hydrothermal effects','gravity residual - real data','gravity changes due to magma','gravity changes for model A')
title('La Pietra')



