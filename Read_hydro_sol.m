function [Pu,Tu,dRhou]=Read_hydro_sol(filename,Xu,Yu,iter,steady)
global porosity
%filename='\\ctmgnas\UFGM\currenti\ThermoPoro\HalfSpace\HS_Het\Plot_scalar.HS_Het';
%steady=11;

Sol=load(filename);
a=find(Sol(:,4)==Sol(1,4));
Ltime=length(a);
Ntime=length(Sol)/Ltime;

T0=Sol((steady-1)*Ltime+1:(steady)*Ltime,5);
P0=Sol((steady-1)*Ltime+1:(steady)*Ltime,6);

T0=Sol((steady-1)*Ltime+1:(steady)*Ltime,5);
P0=Sol((steady-1)*Ltime+1:(steady)*Ltime,6);
Rhog0=1000*Sol((steady-1)*Ltime+1:(steady)*Ltime,9);
Rhof0=1000*Sol((steady-1)*Ltime+1:(steady)*Ltime,8);
Sat0=Sol((steady-1)*Ltime+1:(steady)*Ltime,7);
Rho0=porosity*(Rhof0.*Sat0+Rhog0.*(1-Sat0));

X=Sol(1:Ltime,1)*1e3;
Z=Sol(1:Ltime,3)*1e3-1500;

time=Sol((iter-1)*Ltime+1,4)
T=Sol((iter-1)*Ltime+1:(iter)*Ltime,5)-T0;
P=Sol((iter-1)*Ltime+1:(iter)*Ltime,6)-P0;

Rhog=1000*Sol((iter-1)*Ltime+1:(iter)*Ltime,9);
Rhof=1000*Sol((iter-1)*Ltime+1:(iter)*Ltime,8);
Sat=Sol((iter-1)*Ltime+1:(iter)*Ltime,7);
Rho=porosity*(Rhof.*Sat+Rhog.*(1-Sat));
dRho=Rho-Rho0;

b=min(min(X));
a=find(X==b);
X=[zeros(length(a),1);X];
Z=[Z(a);Z]; T=[T(a);T]; P=[P(a);P]; dRho=[dRho(a);dRho];

Tu=griddata(X,Z,T,Xu,Yu);
Pu=griddata(X,Z,P,Xu,Yu);
dRhou=griddata(X,Z,dRho,Xu,Yu);

a=isnan(Tu);
Tu(a)=0;
a=isnan(Pu);
Pu(a)=0;
a=isnan(dRhou);
dRhou(a)=0;


