function [Pu,Tu,dRhou]=Read_TOUGH_sol(filename,Xu,Yu,iter,steady)

Fault1P1=[2541.54 -3000];
Fault1P2=[3000 -200];
Fault2P1=[5696.16 -3000];
Fault2P2=[6500    0];

% global porosity
% steady=10;
% iter=13;
% filename='output/Sol_';

Sol=load([filename num2str(steady)]);
X=Sol(:,1);
Y=Sol(:,2);
Z=Sol(:,3);
%[porosity,RockDensity]=por_and_rho(X,Y,Z,Fault1P1,Fault1P2,Fault2P1,Fault2P2);
[porosity,RockDensity]=por_and_rho(X,Y,Z);
% [porosity,RockDensity]=Por(X,Y,Z);
T0=Sol(:,4);
P0=Sol(:,5);
Sg0=Sol(:,6);
Dw0=Sol(:,7);
Dg0=Sol(:,8);
Rho0=porosity.*(Dw0.*(Sg0-1)+Dg0.*Sg0);

Sol=load([filename num2str(iter)]);
T=Sol(:,4)-T0;
P=Sol(:,5)-P0;
Sg=Sol(:,6);%-Sg0;
Dw=Sol(:,7);%-Dw0;
Dg=Sol(:,8);%-Dg0;
Rho=porosity.*(Dw.*(Sg-1)+Dg.*Sg); %+(1-porosity).*RockDensity;
dRho=Rho-Rho0;

b=min(min(X));
a=find(X==b);
X=[zeros(length(a),1);X];
Z=[Z(a);Z]; T=[T(a);T]; P=[P(a);P];dRho=[dRho(a);dRho];

% % % Zmin=(3*Z(1)-Z(2))/2;
% % % X=[X;X];
% % % Z=[Z;2*Zmin-Z];
% % % T=[T;T];
% % % P=[P;P];
% % % dRho=[dRho;dRho];

Tu=griddata(X,Z,T,Xu,Yu);
Pu=griddata(X,Z,P,Xu,Yu);
dRhou=griddata(X,Z,dRho,Xu,Yu);

a=isnan(Tu);
Tu(a)=0;
a=isnan(Pu);
Pu(a)=0;
a=isnan(dRhou);
dRhou(a)=0;





