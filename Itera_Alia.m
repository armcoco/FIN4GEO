clear all
close all
clc

global miu Kb alphaalpha betabeta;

h=1;
for miu=[2]*1e9
    for Kb=5e9 %[1 5 10 15 20]*1e9
        for alphaalpha=1
            for betabeta=1
                for fault=0:1
                    LinearElasticity
                    
                    B_disl=0.05;
                    dip=1.768191886644777;
                    trasl=4300;
                    s1=200/cos(dip-pi/2);
                    s2=1500/cos(dip-pi/2);
                    [u1,v1]=InclinedDip3(Xusup-trasl,0,s1,s2,dip,ni,B_disl*2,miu);
                    dip=1.793696559123272;
                    trasl=7940;
                    s1=0/cos(dip-pi/2);
                    s2=1500/cos(dip-pi/2);
                    [u2,v2]=InclinedDip3(Xusup-trasl,0,s1,s2,dip,ni,B_disl*2,miu);
                    u_fault=u1+u2;
                    v_fault=v1+v2;
                    usup=usup+fault*u_fault';
                    vsup=vsup+fault*v_fault';
                    
                    Utime(:,h)=usup;
                    Vtime(:,h)=vsup;
                    Matrix{h}=M;
                    %                 legend_str{h}=['\mu = ',num2str(miu/1e9),' GPa'];
                    %                 legend_str{h}=['K = ',num2str(Kb/1e9),' GPa'];
                    legend_str{h}=['fault = ',num2str(fault)];
                    h=h+1
                end
            end
        end
    end
end


%%
FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

figure(1)
subplot(1,2,1)
plot(Xusup,Utime,'lineWidth',1.1)
xlim([0 10000])
xlabel('radial distance [m]')
title('Horizontal displacement [m]')
subplot(1,2,2)
plot(Xusup,Vtime,'lineWidth',1.1)
xlim([0 10000])
title('Vertical displacement [m]')
xlabel('radial distance [m]')
legend(legend_str)
% legend(num2str([1:5]'))
