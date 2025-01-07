clear all
clc

%%%%%for Deformation
Software='tough';
global Software IfPhi2 alphaalpha expmiu

pathSol='Faults_paper_new/ScenarioA/unrest_x1/';

ToughSolFile=[pathSol,'Sol_'];
tt=load([pathSol,'time'])/(365.25*86400);
Ntime=length(tt)-1;

h=1;
steady=1000;

times_sel=15:25;

selected_times=0:1:Ntime;
for iter=0:2 %selected_times
    expmiu=1;
    alphaalpha=1;
    IfPhi2=0;
    LinearElasticity_header;
    LinearElasticity_main;
    Utime(:,h)=usup;
    Vtime(:,h)=vsup;
    U_axis(h)=u_axis;
    V_axis(h)=v_axis;
    U_axis(h)=u_axis;
    V_axis(h)=v_axis;
    dGtime(:,h)=dGsup*1e8;
    h=h+1;
end

return

close all
figure(1)
subplot(1,2,1)
plot(Xusup,[Utime,UtimePS,UtimePSobl,UtimePSpro])
xlim([0 10000])
xlabel('r[m]')
title('U[m]')
subplot(1,2,2)
plot(Xusup,[Vtime,VtimePS,VtimePSobl,VtimePSpro])
xlim([0 10000])
title('V[m]')
xlabel('r[m]')
legend('poroelastic','with spherical source','with prolate source','with oblate source')
% h1=legend(num2str(tt(times_sel)))
% v1 = get(h1,'title');
% set(v1,'string','time[y]');
return

figure(2)
subplot(1,2,1)
plot(tt(1:Ntime+1)',Utime(1,:),'-o')
xlabel('time[y]')
title('U[m]')
%xlim([0 5])
subplot(1,2,2)
plot(tt(1:Ntime+1)',Vtime(1,:),'-o')
title('V[m]')
xlabel('time[y]')
%xlim([0 5])

figure(3)
subplot(1,2,1)
plot(tt(1:Ntime+1)',U_axis,'-o')
xlabel('time[y]')
title('U[m]')
%xlim([0 5])
subplot(1,2,2)
plot(tt(1:Ntime+1)',V_axis,'-o')
title('V[m]')
xlabel('time[y]')
%xlim([0 5])
return
%%%%%For Gravity

global HydroSolFile;
HydroSolFile='\\ctmgnas\UFGM\currenti\ThermoPoro\CampiFlegrei3\Plot_scalar.CF3Rin1_5km';
Sol=load(HydroSolFile);
a=find(Sol(:,4)==Sol(1,4));
Ntime=length(Sol)/length(a);

h=1;
for iter=12:Ntime
    LinearElasticity
    dGtime(:,h)=dGsupT*1e8;
    h=h+1
end

figure(1)
plot(traccia(:,1),dGtime)
xlim([0 3000])
xlabel('r[m]')
ylabel('\Deltag[\muGal]')

figure(2)
plot([0.1:0.1:3]',dGtime(2,:))
xlabel('time[y]')
ylabel('\Deltag[\muGal]')
xlim([0 3])

