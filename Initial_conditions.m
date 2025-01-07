clear all
close all

steady=1000;
% unrest=10;
% endsim=48;

load Faults_paper_ScenarioA_x1.mat

Sol=load(['Sol_',num2str(steady)]);
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


subplot(1,2,1)
[c,h]=contourf(XX,ZZ,PP/1e6);
clabel(c,h);
% caxis([-0.3 2])
% xlim([0 1000])
colorbar
title('Initial condition: pore pressure [MPa]','FontSize',16)

subplot(1,2,2)
[c,h]=contourf(XX,ZZ,TT);
clabel(c,h);
% caxis([-10 60])
% xlim([0 1000])
colorbar
title('Initial condition: temperature [K]','FontSize',16)
