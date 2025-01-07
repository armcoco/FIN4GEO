clear all
clc
close all

global Software IfPhi2 alphaalpha expmiu c_parameterx c_parameterz

%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathSol='tough2_sol_files/'; %where tough outputs files are located
steady=1000; %e.g. if steady=1000 then Sol1000 is the reference state (steady state) to which compute the pressure and temperature differences at each time step
n=416; % number of grid points per grid direction
c_parameterx=1e4; %reference scale in the x direction (there will be an adequate resolution for 0<x<c_parameterx)
c_parameterz=3e3; %reference scale in the z direction (there will be an adequate resolution for |z|<c_parameterz)
% DEFINE macroscopic streaming current coupling coefficient in L_SP__.m
% DEFINE electrical conductivity of the porous media in sig__.m
% DEFINE RIGIDITY IN G__.m (MIND THE DOUBLE UNDERSCORE)
% DEFINE POISSON RATIO IN nu__.m (MIND THE DOUBLE UNDERSCORE)
% DEFINE TOPOGRAPHY IN phi1__.m (MIND THE DOUBLE UNDERSCORE)
% DEFINE POROSITY AND DENISTY IN por_and_rho.m 
% DEFINE selected_index_times later in this file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Software='tough';

ToughSolFile=[pathSol,'Sol_'];
tt=load([pathSol,'time'])/(365.25*86400);
Ntime=length(tt)-1;

h=1;

selected_index_times=0:1:Ntime;
selected_index_times=[20 48]; %[1,2,21,41,46,49]-1;
for iter=selected_index_times
    disp(['time: ',num2str(tt(iter+1))])
    expmiu=1;
    alphaalpha=1;
    IfPhi2=0;
    LinearElasticity_header;
    LinearElasticity_main;
    Utime(:,h)=[u_axis;usup];
    Vtime(:,h)=[v_axis;vsup];
    U_axis(h)=u_axis;
    V_axis(h)=v_axis;
    U_axis(h)=u_axis;
    V_axis(h)=v_axis;
    dGtime(:,h)=[dG_axis;dGsup]*1e8;
    dGtimeDefo(:,h)=[dG_axisDefo;dGsupDefo]*1e8;
    sptime(:,h)=[sp_axis;spsup];
    h=h+1;
end

xtopo=[0,Xusup];
ztopo=[2440,Yusup];
xgrid=Xu;
zgrid=Yu;

figure(1)
hold on
plot(xgrid,zgrid,'.b')
xlabel('x [m]','FontSize',16)
ylabel('z [m]','FontSize',16)
xlim([0 1e4])
ylim([-1e4 1e4])
contour(Xu,Yu,Phi1,[0 0],'r','LineWidth',2)
title('Grid points and topography (ZOOM IN AND OUT TO SEE MORE OR LESS POINTS)','FontSize',12)

figure(2)
plot(xtopo,Utime)
xlabel('x [m]','FontSize',16)
ylabel('\Delta U [m]','FontSize',16)
title('Horizontal ground displacement at different times')
xlim([0 1e4])

figure(3)
plot(xtopo,Vtime)
xlabel('x [m]','FontSize',16)
ylabel('\Delta V [m]','FontSize',16)
title('Vertical ground displacement at different times')
%xlim([0 5000])
xlim([0 1e4])


figure(4)
plot(xtopo,dGtime)
xlabel('x [m]','FontSize',16)
ylabel('\Deltag [\muGal]','FontSize',16)
title('gravity changes (on the ground) at different times')
xlim([0 1e4])

figure(5)
plot(xtopo,sptime*1000)
xlabel('x [m]','FontSize',16)
ylabel('Self-potential[mV]','FontSize',16)
title('self potential (on the ground) at different times')
xlim([0 1e4])

figure(6)
plot(xtopo,ztopo)
xlabel('x [m]','FontSize',16)
ylabel('z [m]','FontSize',16)
title('topography')
xlim([0 1e4])

figure(7)
plot(xtopo,dGtimeDefo)
xlabel('x [m]','FontSize',16)
ylabel('\Deltag [\muGal]','FontSize',16)
title('gravity changes (on the ground) including defo effect')
xlim([0 1e4])

%Plot 0d, 1d, 10d, 100d, 350d, 500d
sptime2=sptime(:,[1,2,21,41,46,49]);
Utime2=Utime(:,[1,2,21,41,46,49]);
Vtime2=Vtime(:,[1,2,21,41,46,49]);
dGtime2=dGtime(:,[1,2,21,41,46,49]);


figure(7)
plot(xtopo,Utime2)
xlabel('Radial distance [m]','FontSize',16)
ylabel('\Delta U [m]','FontSize',16)
legend('0d','1d','10d','100d','350d','500d')
title('Horizontal ground displacement at different times')
xlim([0 5500])

figure(8)
plot(xtopo,Vtime2)
xlabel('Radial distance [m]','FontSize',16)
ylabel('\Delta V [m]','FontSize',16)
title('Vertical ground displacement at different times')
legend('0d','1d','10d','100d','350d','500d')
%xlim([0 5000])
xlim([0 5500])


figure(9)
plot(xtopo,dGtime2)
xlabel('Radial distance [m]','FontSize',16)
ylabel('\Deltag [\muGal]','FontSize',16)
title('Gravity changes (on the ground) at different times')
legend('0d','1d','10d','100d','350d','500d')
xlim([0 5500])


figure(10)
plot(xtopo,sptime2*1000)
xlabel('Radial distance [m]','FontSize',16)
ylabel('Self-potential[mV]','FontSize',16)
title('Self-potential (on the ground) at different times')
legend('0d','1d','10d','100d','350d','500d')
xlim([0 5500])
