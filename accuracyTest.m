clear all
close all
clc

global n c_parameterx c_parameterz academic_test;

c_parameterx=1e4; %reference scale in the x direction (there will be an adequate resolution for 0<x<c_parameterx)
c_parameterz=3e3; %reference scale in the z direction (there will be an adequate resolution for |z|<c_parameterz)

academic_test = 1;
Software='tough';

N = ceil([100]*1.2.^[0:9]);
N = ceil([64]*1.2.^[0:6]);
N = [50 100];
% N = ceil([32]*1.2.^[0:6]);
% N=[100];
% N=[50:25:100];
% N=[50:5:110];
% N=[50:70];

h=1;
for n=N
    n
%     mainMG_rho;
%     RHO(:,kk)=Rho;
%     RHO_MEAN(:,kk)=Rho_mean;
%     RHO_MEAN_10(:,kk)=Rho_mean_10;
%     mainPS;
    expmiu=1;
    alphaalpha=1;
    IfPhi2=0;
    LinearElasticity_header;
    LinearElasticity_main;
%     Utime(:,h)=usup;
%     Vtime(:,h)=vsup;
%     U_axis(h)=u_axis;
%     V_axis(h)=v_axis;
%     U_axis(h)=u_axis;
%     V_axis(h)=v_axis;
%     dGtime(:,h)=dGsup*1e8;
%     sptime(:,h)=spsup;
    h=h+1;
end

% fid=fopen('dati_flower_1e0_1e6.txt','w');
% fprintf(fid,'%d \n',length(N));
% fprintf(fid,'%d \n',N);
% fprintf(fid,'%.16f \n',LL1S);
% fprintf(fid,'%.16f \n',LLiS);
% fprintf(fid,'%.16f \n',LL1G);
% fprintf(fid,'%.16f \n',LLiG);
% fclose(fid);
