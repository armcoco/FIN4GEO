clear all
close all
clc

load Serapeo.txt
% load Grav_paper_new_ScenarioA_x0p5.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x0p5.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1.mat
% load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17_hybrid.mat
load Grav_paper_new_extended_release_2km_ScenarioA_x1p15_molar_ratio_0p17_linear_decay.mat
dG_axis=dGtime(1,:);
dG_axis=[0,dG_axis];
% dG_axis=max(0,dG_axis);
delay=1*1.75;
ttt=[0;tt];
Ser=interp1(Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,ttt);

plot(ttt,-dG_axis)
hold on
xlim([0 20])
plot(ttt,Ser)
plot(Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,'*')
df=Ser+dG_axis';
plot(ttt,df,'k')
plot(ttt,470*(1-exp(-0.2*ttt))','r')

plot(ttt,470*(1-exp(-0.2*ttt))-dG_axis(:),'r')


