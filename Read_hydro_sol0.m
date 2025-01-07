filename='\\ctmgnas\UFGM\currenti\ThermoPoro\HalfSpace\HS_Het\Plot_scalar.HS_Het';
Sol=load(filename);
a=find(Sol(:,4)==Sol(1,4));
Ltime=length(a);
Ntime=length(Sol)/Ltime;
porosity=0.2;

steady=11;
T0=Sol((steady-1)*Ltime+1:(steady)*Ltime,5);
P0=Sol((steady-1)*Ltime+1:(steady)*Ltime,6);

Rhog0=1000*Sol((steady-1)*Ltime+1:(steady)*Ltime,9);
Rhof0=1000*Sol((sxlteady-1)*Ltime+1:(steady)*Ltime,8);
Sat0=Sol((steady-1)*Ltime+1:(steady)*Ltime,7);
Rho0=porosity*(Rhof0.*Sat0+Rhog0.*(1-Sat0));

X=Sol(1:Ltime,1)*1e3;
Z=Sol(1:Ltime,3)*1e3-1500;


for iter=1:Ntime
    T(:,iter)=Sol((iter-1)*Ltime+1:(iter)*Ltime,5);
    time(:,iter)=Sol((iter-1)*Ltime+1:(iter)*Ltime,4);
end


for iter=12:53
time(iter-steady)=Sol((iter-1)*Ltime+1,4);
T(:,iter-steady)=Sol((iter-1)*Ltime+1:(iter)*Ltime,5)-T0;
P(:,iter-steady)=Sol((iter-1)*Ltime+1:(iter)*Ltime,6)-P0;
% dRho(:,iter-steady)=(1000*Sol((iter-1)*Ltime+1:(iter)*Ltime,8)-dRho0)*0.2;
% dRho1(:,iter-steady)=0.2*1000*(Sol((iter-1)*Ltime+1:(iter)*Ltime,8)-Sol((iter-1)*Ltime+1:(iter)*Ltime,9)).*(Sol((iter-1)*Ltime+1:(iter)*Ltime,7)-1);

Rhog(:,iter-steady)=1000*Sol((iter-1)*Ltime+1:(iter)*Ltime,9);
Rhof(:,iter-steady)=1000*Sol((iter-1)*Ltime+1:(iter)*Ltime,8);
Sat(:,iter-steady)=Sol((iter-1)*Ltime+1:(iter)*Ltime,7);

Rho(:,iter-steady)=porosity*(Rhof(:,iter-steady).*Sat(:,iter-steady)+Rhog(:,iter-steady).*(1-Sat(:,iter-steady)));
dRho(:,iter-steady)=Rho(:,iter-steady)-Rho0;

MassFraction(:,iter-steady)=Sat(:,iter-steady).*Rhof(:,iter-steady)./(Sat(:,iter-steady).*Rhof(:,iter-steady)+(1-Sat(:,iter-steady)).*Rhog(:,iter-steady));
end

%dRho=phi*(drho_f-drho_g)*(S-1)
%dRho=porosity*1000*(Sol((iter-1)*Ltime+1:(iter)*Ltime,8)-Sol((iter-1)*Ltime+1:(iter)*Ltime,9)).*(Sol((iter-1)*Ltime+1:(iter)*Ltime,7)-1);

Tmin=min(min(T));
Tmax=max(max(T));

Pmin=min(min(P));
Pmax=max(max(P));

for j=1:42
surf(reshape(X,60,30),reshape(Z,60,30),reshape(T(:,j),60,30))
%[c,h]=contour(reshape(X,60,30),reshape(Z,60,30),reshape(T(:,j),60,30),[Tmin:floor(Tmax-Tmin)/10:Tmax]);
clabel(c,h)
title(num2str(tt(j)-10000))
view(2)
shading interp
caxis([Tmin Tmax])
colorbar
axis equal
xlim([0 2000])
saveas(gcf,['Temp\' num2str(j)],'jpg')

surf(reshape(X,60,30),reshape(Z,60,30),reshape(P(:,j),60,30))
title(num2str(tt(j)-10000))
view(2)
caxis([Pmin Pmax])
colorbar
axis equal
xlim([0 2000])
saveas(gcf,['Pres\' num2str(j)],'jpg')
end


Tu=griddata(X,Z,T,Xu,Yu);
Pu=griddata(X,Z,P,Xu,Yu);
dRhou=griddata(X,Z,dRho,Xu,Yu);

a=isnan(Tu);
Tu(a)=0;
a=isnan(Pu);
Pu(a)=0;
a=isnan(dRhou);
dRhou(a)=0;

return
%%% ControlloSteady
for iter=1:10
times(iter)=Sol((iter-1)*Ltime+1,4);
Ts(:,iter)=Sol((iter-1)*Ltime+1:(iter)*Ltime,5);
Ps(:,iter)=Sol((iter-1)*Ltime+1:(iter)*Ltime,6);
end



