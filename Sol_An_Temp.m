function [Ux,Uz]=Sol_An_Temp(x,y,Ra,Rb,Ta,miu)
global alpha
lambda=miu;
Par=alpha*(3*lambda+2*miu)/(lambda+2*miu);
r=sqrt(x.^2+y.^2);


%Linear Temperature Function in a shell Radius and Rb
At=Ta/(Ra-Rb);
Bt=-Ta*Rb/(Ra-Rb);
u=Par./r.^2.*(IntTemp(r,At,Bt)-IntTemp(Ra,At,Bt)).*(r<=Rb)+Par./r.^2.*(IntTemp(Rb,At,Bt)-IntTemp(Ra,At,Bt)).*(r>Rb);

%Parabolic Temperature Function
u=Par./r.^2.*Ta*1e4.*(r-Ra);


Ux=u.*x./r;
Uz=u.*y./r;

function F=IntTemp(x,At,Bt)
F=At*x.^4/4+Bt*x.^3/3;
return