clear all
close all

FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);


steady=0;
% unrest=10;
% endsim=48;

filename='coupling_with_MUFITS/SCENARIO1/SCENARIO1';

filename_grid=[filename,'.GRID.SUM'];
filename_out_ref=[filename,'.',num2str(steady,'%.4d'),'.SUM'];

[startRow, endRow]=read_from_Mufits(filename_grid);
dataArray=read_data_from_Mufits(filename_grid,startRow,endRow,5);
cellID_grid=dataArray{:,1};
X=dataArray{:,2};
Y=dataArray{:,3};
Z=-dataArray{:,4};
porosity=dataArray{:,5};

[startRow, endRow]=read_from_Mufits(filename_out_ref);
dataArray=read_data_from_Mufits(filename_out_ref,startRow,endRow,4);
cellID_out0=dataArray{:,1};
P0=1e6*dataArray{:,2};
T0=dataArray{:,3};
Rho0=dataArray{:,4};
Rho0=porosity.*Rho0;

a=find(Z==Z(1));
l=length(a);
m=length(Z)/l;
TT0=reshape(T0,l,m);
PP0=reshape(P0,l,m);
RRho0=reshape(Rho0,l,m);

% b=min(min(X));
% a=find(X==b);
% X=[zeros(length(a),1);X];
% Z=[Z(a);Z]; Y=[Y(a);Y]; dT=[dT(a);dT]; dP=[dP(a);dP];dRho=[dRho(a);dRho];

time=loadMufitsTimes(filename);
for iter=0:8 %[1,11,39,59] %0:length(time)-1
    filename_out=[filename,'.',num2str(iter,'%.4d'),'.SUM'];
    time_=time(iter+1)/365.25;
    clf
    
    dataArray=read_data_from_Mufits(filename_out,startRow,endRow,4);
    cellID_out=dataArray{:,1};
    P=1e6*dataArray{:,2};
    T=dataArray{:,3};
    Rho=dataArray{:,4};
    Rho=porosity.*Rho;

    if max(abs(cellID_grid-cellID_out))
        disp('WARNING!!! cellID_grid and cellID_out are differents!!!')
        qqq
        return;
    end

    a=find(Z==Z(1));
    l=length(a);
    m=length(Z)/l;
    XX=reshape(X,l,m);
    ZZ=reshape(Z,l,m);
    TT=reshape(T,l,m);
    PP=reshape(P,l,m);
    RRho=reshape(Rho,l,m);
    
    subplot(3,1,1)
    [c,h]=contourf(XX,ZZ,(PP-PP0)/1e6);
%     clabel(c,h);
%     caxis([-0.5 30])
% %     caxis([0 4])
    %     xlim([0 1000])
    colorbar
    title(['pore pressure changes P-P_0 [MPa],   time: ',num2str(round(time_*100)/100),' years'],'FontSize',16)
    
    subplot(3,1,2)
    [c,h]=contourf(XX,ZZ,TT-TT0);
%     clabel(c,h);
%     caxis([0 360])
% %     caxis([-15 100])
% %     caxis([-15 170])
    %     xlim([0 1000])
    colorbar
    title('temperature changes T-T_0 [C]','FontSize',16)
    hold on
%     plottaGradT;
    
    subplot(3,1,3)
    [c,h]=contourf(XX,ZZ,-(RRho-RRho0));
%     [c,h]=contourf(XX,ZZ,SSGG-SSGG0);
%     clabel(c,h);
%     caxis([0 360])
% %     caxis([-15 200])
% %     caxis([-15 320])
% %     caxis([0 1])
    %     xlim([0 1000])
    colorbar
%     title('density changes -(\rho-\rho_0) [Kg/m^3]','FontSize',16)
    title('gas saturation changes S_g-S_{g0}','FontSize',16)
    
    drawnow;
    disp(['time: ',num2str((time(iter+1)/365.25)),' years'])
    pause;
end
