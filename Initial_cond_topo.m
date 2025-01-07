clear all
close all

steady=0;
% unrest=10;
% endsim=48;

% load Faults_paper_ScenarioA_x1.mat

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


[c,h]=contourf(XX,ZZ,TT);
clabel(c,h);
% caxis([-10 60])
% xlim([0 1000])
colorbar
title('Initial condition: temperature [C]','FontSize',16)

