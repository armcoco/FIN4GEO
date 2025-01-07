clear all
close all

FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


steady=1000;
% unrest=10;
% endsim=48;

Sol=load(['Alia/HighInj_increasedRate/Sol_',num2str(steady)]);
X=Sol(:,1);
Y=Sol(:,2);
Z=Sol(:,3);
a=find(Z==Z(1));
l=length(a);
m=length(Z)/l;
XX=reshape(X,l,m);
ZZ=reshape(Z,l,m);
porosity=Por(X,Y,Z);
Sg0=Sol(:,6);
Dw0=Sol(:,7);
Dg0=Sol(:,8);
Rho0=porosity.*(Dw0.*(1-Sg0)+Dg0.*Sg0);
RRho0=reshape(Rho0,l,m);
SSg0=reshape(Sg0,l,m);

time=load('Alia/HighInj/time');
for iter=0:length(time)-1
    clf
    Sol=load(['Alia/HighInj/Sol_',num2str(iter)]);
    Sg=Sol(:,6);
    Dw=Sol(:,7);
    Dg=Sol(:,8);
    Rho=porosity.*(Dw.*(1-Sg)+Dg.*Sg);
    
    RRho=reshape(Rho,l,m);
    SSg=reshape(Sg,l,m);
    
    [c,h]=contourf(XX,ZZ,RRho-RRho0);
%         [c,h]=contourf(XX,ZZ,SSg);
    clabel(c,h);
    % caxis([-0.3 2])
    xlim([0 10000])
    colorbar
%     title('pore pressure [MPa]','FontSize',16)
    
    drawnow;
    disp(['time: ',num2str(time(iter+1)/86400/365.25)])
    pause;
end
