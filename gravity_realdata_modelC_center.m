%%% model A

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
% load Grav_paper_new_extended_release_200m_ScenarioA_x1p15_molar_ratio_0p17_Inj2000m_hybrid.mat
load Grav_paper_new_extended_release_200m_ScenarioA_x1p15_molar_ratio_0p17_Inj2000m_hybrid_2.mat
dG_axis=dGtime(20,:);
dG_axis=[0,dG_axis];
% dG_axis=max(0,dG_axis);
delay=1*1.75;
delay=1*2.10;
ttt=[0;tt]+1982;
Ser=interp1(1982+Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,ttt);

plot(ttt,-dG_axis)
hold on
xlim(1982+[0 15])
% plot(ttt,Ser)
plot(1982+Serapeo(:,1)/365.25-delay,Serapeo(:,2)*1e3,'-*')
df=Ser+dG_axis';
% plot(ttt,df,'y')
defo_magma=320*(1-exp(-0.26*(ttt-1982)));
plot(ttt,defo_magma','k')
plot(ttt,defo_magma-dG_axis(:),'r')

legend('gravity changes due to hydrothermal effects','gravity residual - real data','gravity changes due to magma','gravity changes for model A')
title('Serapeo')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load LaPietra.txt
figure
factor=3.4;
dG_axis=factor*dGtime(40,:);
dG_axis=[0,dG_axis];
% dG_axis=max(0,dG_axis);
delay=1*1.75;
delay=1*2.10;
ttt=[0;tt]+1982;
Solf=interp1(1982+LaPietra(:,1)/365.25-delay,LaPietra(:,2)*1e3,ttt);

plot(ttt,-dG_axis)
hold on
xlim(1982+[0 15])
% plot(ttt,Ser)
plot(1982+LaPietra(:,1)/365.25-delay,LaPietra(:,2)*1e3,'-*')
df=Solf+dG_axis';
% plot(ttt,df,'y')
% defo_magma=360*(1-exp(-0.16*ttt));
defo_magma=factor*72*(1-exp(-0.22*(ttt-1982)));
plot(ttt,defo_magma','k')
plot(ttt,defo_magma-dG_axis(:),'r')

legend('gravity changes due to hydrothermal effects','gravity residual - real data','gravity changes due to magma','gravity changes for model A')
title('La Pietra')

figure 
magma_acc=(1.80+0.36)*(1-exp(-0.24*(ttt-1982)));
plot(ttt,magma_acc)
xlim([0 15]+1982)
ylabel('magma volme increase [km^3]')

