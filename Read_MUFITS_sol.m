function [Pu,Tu,dRhou]=Read_MUFITS_sol(filename,Xu,Yu,iter,steady)

filename_grid=[filename,'.GRID.SUM'];
filename_out=[filename,'.',num2str(iter,'%.4d'),'.SUM'];
filename_out_ref=[filename,'.',num2str(steady,'%.4d'),'.SUM'];

[startRow, endRow]=read_from_Mufits(filename_grid);
dataArray=read_data_from_Mufits(filename_grid,startRow,endRow,5);
cellID_grid=dataArray{:,1};
X=dataArray{:,2};
Y=dataArray{:,3};
Z=-dataArray{:,4};
porosity=dataArray{:,5};

[startRow, endRow]=read_from_Mufits(filename_out);
dataArray=read_data_from_Mufits(filename_out,startRow,endRow,4);
cellID_out=dataArray{:,1};
P=1e6*dataArray{:,2};
T=dataArray{:,3};
Rho=dataArray{:,4};
dataArray=read_data_from_Mufits(filename_out_ref,startRow,endRow,4);
cellID_out0=dataArray{:,1};
P0=1e6*dataArray{:,2};
T0=dataArray{:,3};
Rho0=dataArray{:,4};

if max(abs(cellID_grid-cellID_out))
    disp('WARNING!!! cellID_grid and cellID_out are differents!!!')
    qqq
    return;
end

Rho0=porosity.*Rho0;
Rho=porosity.*Rho;

dP=P-P0;
dT=T-T0;
dRho=Rho-Rho0;

b=min(min(X));
a=find(X==b);
X=[zeros(length(a),1);X];
Z=[Z(a);Z]; Y=[Y(a);Y]; dT=[dT(a);dT]; dP=[dP(a);dP];dRho=[dRho(a);dRho];

% % % Zmin=(3*Z(1)-Z(2))/2;
% % % X=[X;X];
% % % Z=[Z;2*Zmin-Z];
% % % T=[T;T];
% % % P=[P;P];
% % % dRho=[dRho;dRho];

Tu=griddata(X,Z,dT,Xu,Yu);
Pu=griddata(X,Z,dP,Xu,Yu);
dRhou=griddata(X,Z,dRho,Xu,Yu);

a=isnan(Tu);
Tu(a)=0;
a=isnan(Pu);
Pu(a)=0;
a=isnan(dRhou);
dRhou(a)=0;

