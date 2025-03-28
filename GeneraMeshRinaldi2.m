clear
close all

%Distribuzione serie geometrica
L=10000-320; %lunghezza massima
q=1.075; %fatttore
N=49;
l=L*(1-q)/(1-q^N);
sx=l*q.^[0:1:N-1]; 
sx=[20*ones(16,1);sx'];

%Distribuzione simile a flow.inp con faglie
%sx=[20*ones(40,1);200*ones(15,1);20*ones(40,1);200*ones(15,1);20*ones(109,1);200*ones(2,1)];

Tratti=[0;cumsum(sx)];

for i=1:length(Tratti)-1
    X(i)=(Tratti(i)+Tratti(i+1))/2;
end

stepV=20;
layer=-1490.0;
top=-10.0;
Row=(-layer+top)/stepV+2;
Col=length(X);
D=1e-9;
LL=[layer:stepV:top top+1]';


for i=1:length(X)
    Volume(i,1)=pi*(Tratti(i+1)^2-Tratti(i)^2)*stepV;
end

Xcenter=[];
Zcenter=[];
for i=1:Row  
    Xcenter=[Xcenter;X'];
    Zcenter=[Zcenter; LL(i)*ones(size(X'))];   
end
Ycenter=0.5*ones(size(Xcenter));
figure(1)
plot(Xcenter,Zcenter,'o')
save xcenter.txt Xcenter -ascii -double
save ycenter.txt Ycenter -ascii -double
save zcenter.txt Zcenter -ascii -double
clear Xcenter Ycenter Zcenter

fid = fopen('flow.inp','w');

fprintf(fid,'TOUGH2 Analysis\n');
fprintf(fid,'ROCKS----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'Rock1    3    2000.0       0.2   1.0E-14   1.0E-14   1.0E-14      2.80    1000.0\n');
fprintf(fid,'       0.0       0.0      2.80       0.0       0.0\n');
fprintf(fid,'    1            0.2       0.1       0.9       0.7\n');
fprintf(fid,'    1            0.0       0.0       1.0\n');
fprintf(fid,'\n');
fprintf(fid,'MULTI----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'    2    3    2    6\n');
fprintf(fid,'START----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'PARAM----1-MOP* 123456789012345678901234----*----5----*----6----*----7----*----8\n');
fprintf(fid,' 8 29999    1000000000000002  03 000   0                                        \n');
fprintf(fid,'       0.0   1.27E11     100.0   3.15E10                9.81       4.0       1.0\n');
fprintf(fid,'    1.0E-5       1.0                 1.0       1.0          \n');
fprintf(fid,'\n');
fprintf(fid,'SOLVR----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'5  Z1   O0       0.1    1.0E-6\n');

%%%%%  FILE for TIMES  %%%%%%%%%%%%%%%%%%
times=[[0:2000:4000]';[4000.1:0.5:4010]']*86400*365';
l=length(times);
fprintf(fid,'TIMES----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'%5d%5d\n',l,l); 
h=1;
i=1;
while i<=l
    while h<=8 & i<=l
        fprintf(fid,'%#10.5G',times(i));
        h=h+1;
        i=i+1;
    end
    fprintf(fid,'\n');
    h=1;
end

%%%%%  FILE for ELEMENT  %%%%%%%%%%%%%%%%%%
h=1;
fprintf(fid,'ELEME----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
for j=1:Row    
    for i=1:Col
        if j==Row;
            fprintf(fid,'%5d          Rock1%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G\n',(Col)*(h-1)+i, 1e59, 0.0, 1.0, X(i),0.5, LL(j));
        else
            fprintf(fid,'%5d          Rock1%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G\n',(Col)*(h-1)+i, Volume(i), 0.0, 1.0, X(i),0.5, LL(j));
        end            
    end
    h=h+1;
end
fprintf(fid,'\n');


%%%%%  FILE for CONNETTIVITY  %%%%%%%%%%%%%%%%%%
AreaLat=2*pi*Tratti(2:end)*stepV;
%AreaLatTop=2*pi*Tratti(2:end)*D;

for i=1:length(X)
    AreaVer(i)=pi*(Tratti(i+1)^2-Tratti(i)^2);
end

fprintf(fid,'CONNE----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
for i=1:Row    
    for j=1:Col-1
        fprintf(fid,'%5d%5d%20d%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G\n',(i-1)*(Col)+j,(i-1)*(Col)+j+1, 1, Tratti(j+1)-X(j), X(j+1)-Tratti(j+1), AreaLat(j),0, 0);                   
    end
    h=h+1;
end

for i=1:Row-2   
    for j=1:Col
        fprintf(fid,'%5d%5d%20d%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G\n',(i-1)*(Col)+j,(i)*(Col)+j, 1, stepV/2, stepV/2, AreaVer(j),-1, 0);                   
    end
    h=h+1;
end

i=Row-1;
for j=1:Col
    fprintf(fid,'%5d%5d%20d%#10.4G%#10.4G%#10.4G%#10.4G%#10.4G\n',(i-1)*(Col)+j,(i)*(Col)+j, 1, stepV/2, D, AreaVer(j),-1, 0);                   
end


fprintf(fid,'\n');

%%%%%  FILE for INITIAL CONDITION  %%%%%%%%%%%%%%%%%%
fprintf(fid,'INCON----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
GradP=1000*9.81; %in Pa/m
Ptop=101325; %in Pa
GradT=130/1000;
Ttop=20;

h=1;
for i=1:Row
    for j=1:Col
        fprintf(fid,'%5d\n',h);  
        P_CO2=0.1;         if (i==Row) P_CO2=0.0; end
        %fprintf(fid,'%#20.14G%#20.14G%#20.14G\n',-LL(i)*GradP+Ptop,-LL(i)*GradT+Ttop, P_CO2);
        fprintf(fid,'%#20.14G%#20.14G%#20.14G\n',-LL(i)*GradP+Ptop,-LL(i)*GradT+Ttop, Ptop*0.0004);

        h=h+1;
    end   
end
fprintf(fid,'\n');

%%%%%  FILE for GENER  %%%%%%%%%%%%%%%%%%
fprintf(fid,'GENER----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');

%Rinaldi
AreaTot=pi*160^2;
H2O_0=2400*1e3/86400/AreaTot;
H2O_1=6100*1e3/86400/AreaTot;
CO2_0=1000*1e3/86400/AreaTot;
CO2_1=6000*1e3/86400/AreaTot;
Heat=31415926.5/(pi*Tratti(end)^2);
ent_water=2800000;
ent_CO2=290000; %797482;
% ent_water=2980000.0; %usato da Alia
% ent_CO2=115000; %usato da Alia

%%%%%    WATER & CO2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1;
for j=1:8
    fprintf(fid,'%5d%5d                   4     COM11\n',j,h);         
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',0,86400*365*4000.01,86400*365*4001.6,86400*365*4012); %times
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',H2O_0*AreaVer(j),H2O_1*AreaVer(j),H2O_0*AreaVer(j),H2O_0*AreaVer(j)); %rate
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',ent_water,ent_water,ent_water,ent_water); %enthalpy
    h=h+1;
    fprintf(fid,'%5d%5d                   4     COM21\n',j,h);         
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',0,86400*365*4000.01,86400*365*4001.6,86400*365*4012); %times
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',CO2_0*AreaVer(j),CO2_1*AreaVer(j),CO2_0*AreaVer(j),CO2_0*AreaVer(j)); %rate
    fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',ent_CO2,ent_CO2,ent_CO2,ent_CO2); %enthalpy
    h=h+1;                       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%    WATER & CO2 --- enthalpy per m^2   %%%%%%%%%%%%%
% % % h=1;
% % % for j=1:5
% % %     fprintf(fid,'%5d%5d                   4     COM11\n',j,h);         
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',0,86400*365*4000.01,86400*365*4001.6,86400*365*4012); %times
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',H2O_0*AreaVer(j),H2O_1*AreaVer(j),H2O_0*AreaVer(j),H2O_0*AreaVer(j)); %rate
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',ent_water/AreaTot*AreaVer(j),ent_water/AreaTot*AreaVer(j),ent_water/AreaTot*AreaVer(j),ent_water/AreaTot*AreaVer(j)); %enthalpy
% % %     h=h+1;
% % %     fprintf(fid,'%5d%5d                   4     COM21\n',j,h);         
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',0,86400*365*4000.01,86400*365*4001.6,86400*365*4012); %times
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',CO2_0*AreaVer(j),CO2_1*AreaVer(j),CO2_0*AreaVer(j),CO2_0*AreaVer(j)); %rate
% % %     fprintf(fid,'%#14.7G%#14.7G%#14.7G%#14.7G\n',ent_CO2/AreaTot*AreaVer(j),ent_CO2/AreaTot*AreaVer(j),ent_CO2/AreaTot*AreaVer(j),ent_CO2/AreaTot*AreaVer(j)); %enthalpy
% % %     h=h+1;                       
% % % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%     ONLY WATER    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h=1;
% for j=1:8
%     fprintf(fid,'%5d%5d                   3     COM11\n',j,h);         
%     fprintf(fid,'%#14.7G%#14.7G%#14.7G\n',0,86400*365*10000.01,86400*365*10004); %times
%     fprintf(fid,'%#14.7G%#14.7G%#14.7G\n',H2O_0*AreaVer(j),H2O_1*AreaVer(j),H2O_1*AreaVer(j)); %rate
%     fprintf(fid,'%#14.7G%#14.7G%#14.7G\n',2980000.0,2980000.0,2980000.0); %enthalpy
%     h=h+1;                       
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% HEAT INPUT at the BOTTOM %%%%%%%%%%% 
% for j=1:Col
%     fprintf(fid,'%5d%5d                   0     HEAT %#10.4G\n',j,h,Heat*AreaVer(j));         
%     h=h+1;                      
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'\n');
fprintf(fid,'NOVER----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
fprintf(fid,'ENDCY\n');

fclose(fid);

% %Esempio Flow.inp convergente
% AreaTot=pi*200^2;
% H2O_0=13.94/AreaTot; %kg/s
% CO2_0=5.79/AreaTot;
% Heat=31415926.5/(pi*Tratti(end)^2);
% 
% h=1;
% for j=1:8
%     fprintf(fid,'%5d%5d                   0     COM1 %#10.4G%#10.4G\n',j,h,H2O_0*AreaVer(j),2980000.0);         
%     h=h+1;
%     fprintf(fid,'%5d%5d                   0     COM2 %#10.4G%#10.4G\n',j,h,CO2_0*AreaVer(j),115000.0);         
%     h=h+1;                       
% end
% 
% for j=1:Col
%     fprintf(fid,'%5d%5d                   0     HEAT %#10.4G\n',j,h,Heat*AreaVer(j));         
%     h=h+1;                      
% end
% fprintf(fid,'\n');
% fprintf(fid,'NOVER----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8\n');
% fprintf(fid,'ENDCY\n');
% 
% fclose(fid);

