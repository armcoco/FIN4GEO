clear all
% close all
clc

%%%%%for Deformation
Software='tough'; % tough or hydro
global HydroSolFile Software IfPhi2 Rx Ry alphaalpha PPP GGG expmiu

if strcmp(Software,'tough')
    %     pathSol='Faults_Alia/Fault_inj_massflux_per_unit_area_constant/unrest_x1/';
    %     pathSol='Faults_3Km/ScenarioA/unrest_x1/';
    %     pathSol='Alia_flow.out_original/HighInj/';
    %     pathSol='Faults_Alia_same_grid_as_Alia/unrest/';
    %     pathSol='../MultigridAxiDefoGrav 18nov14'/Alia_flow.inp_original/HighInj_withoutINCON_on_flow.inp/INCON_from_steadystate for Low_Fumar_1e-14';
    %     pathSol='Faults/ScenarioC/unrest_x1/';
    %     pathSol='Faults_perm_xi_not_lower/ScenarioB/unrest_x1/';
    %     pathSol='Faults_perm_xi_not_lower_eta_1213/ScenarioB/unrest_x1/';
    %     pathSol='Faults_larger_damage_zone/ScenarioB/unrest_x1/';
    %     pathSol='INCON_from_steadystate for Low_Fumar_1e-14/';
    %     pathSol='INCON_from_SteadyState_Fumar_1e-13/';
    %     pathSol='HighInj/';
    %     pathSol='Alia_HighInj/';
    %     pathSol='Faults_1313_Fumar200/ScenarioC/unrest_x1/';
    %     pathSol='Faults_old/unrest_x1/';
    %     pathSol='Faults_old_first/unrest_x1/';
    %     pathSol='Faults_old_gridFiner/unrest_x1/';
    %     pathSol='Faults_old_only1time/unrest_x1/';
    %     pathSol='Faults_old_vertical_fault/unrest_x1/';
    %     pathSol='another/';
    %     pathSol='Faults_paper/ScenarioA/unrest_x1/';
    %     pathSol='Faults_paper/ScenarioC/unrest_x1_porosity/';
    %     pathSol='Faults_paper_3Km/ScenarioA/unrest_x1/';
%     pathSol='Faults_paper_new/ScenarioA/unrest_x1/';
        pathSol='Faults_paper_new/ScenarioA_2yrs/unrest_x1/';
%         pathSol='Faults_paper_new_TodescoHET2/ScenarioA_2yrs/unrest_x1/';
    ToughSolFile=[pathSol,'Sol_'];
    tt=load([pathSol,'time'])/(365.25*86400);
    Ntime=length(tt)-1;
else
    HydroSolFile='\\ctmgnas\UFGM\currenti\ThermoPoro\HalfSpace\HS_Hom\Plot_scalar.HS_Hom';
    Sol=load(HydroSolFile);
    a=find(Sol(:,4)==Sol(1,4));
    Ntime=length(Sol)/length(a);
    tt=Sol(1:length(a):end,4);
end

% GG0=load('G0sup.txt');

h=1;
steady=1000;

times_sel=15:25; %[12 24 36 48 60];

% radial=load('R.txt');
% load Uzmean.txt
% load Urmean.txt
% load dUzmean.txt
% load dUrmean.txt

PPP=247e6;
% for iter=([2,12,40,60]-1) %[1,3,11,19,39,50,Ntime] %0:Ntime %[1,3,11,19,39] %0:Ntime
V_maxC=[];
V_maxFA=[];
V_maxFB=[];
selected_times=0:5:Ntime;
for iter=selected_times %[11,59]
    PPP=press_time(iter);
%     PPP=press_time(tt(iter+1));
    for expmiu=1 %3.2:0.2:6
        % for iter=0:Ntime
        %     for PPP=6.5e6
        %     Rx=1549;
        %     Ry=Rx/2;
        for alphaalpha=1
%             IfPhi2=0;
%             Software
%             LinearElasticity_header;
%             LinearElasticity_main;
%             %         LinearElasticity;
%             Software
%             if alphaalpha==1
%                 Utime(:,h)=usup;
%                 Vtime(:,h)=vsup;
%                 U_axis(h)=u_axis;
%                 V_axis(h)=v_axis;
%                 U_axis(h)=u_axis;
%                 V_axis(h)=v_axis;
%                 V_maxC(end+1)=v_axis;
%                 V_maxFA(end+1)=max(abs(vsup(abs(Xusup-3000)<1000)));
%                 V_maxFB(end+1)=max(abs(vsup(abs(Xusup-6500)<1000)));
%             else
%                 Utime_onlyP(:,h)=usup;
%                 Vtime_onlyP(:,h)=vsup;
%                 U_axis_onlyP(h)=u_axis;
%                 V_axis_onlyP(h)=v_axis;
%                 U_axis_onlyP(h)=u_axis;
%                 V_axis_onlyP(h)=v_axis;
%             end
                            IfPhi2=1;
%                             Rx=1000;
%                             Ry=Rx;
                            [Rx,Ry]=RxRy_time(tt(iter+1));
                            RhoMagma=RhoMagma_time(tt(iter+1));
                            Software
                            LinearElasticity_header;
                            LinearElasticity_main;
                            Software
%                             Utime(:,h)=usup;
%                             Vtime(:,h)=vsup;
%                             U_axis(h)=u_axis;
%                             V_axis(h)=v_axis;
                            dGtime(:,h)=dGsup*1e8;  
            %     Rx=400;
            %     Ry=800;
            %     Software
            %     LinearElasticity
            %     Software
            %     UtimePSpro(:,h)=usup;
            %     VtimePSpro(:,h)=vsup;
            %     U_axis(h)=u_axis;
            %     V_axis(h)=v_axis;
            %     Rx=1400;
            %     Ry=700;
            %     Software
%                                 LinearElasticity_header;
%                                 LinearElasticity_main;
            %     Software
            %     UtimePSobl(:,h)=usup;
            %     VtimePSobl(:,h)=vsup;
            %     U_axis(h)=u_axis;
            %     V_axis(h)=v_axis;
%                                 dGtime(:,h)=dGsup*1e8;
            %     gradG(:,h)=(dGforGrad*1e8-GG0(2:end-1))./vsup;
        end
        h=h+1
        %     end
    end
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

