clear all
close all
clc

global miu Kb alphaalpha betabeta;

time=load('data/PyTough_output/time_HighInj');
h=1;
for miu=[2]*1e9
    for Kb=5e9 %[1 5 10 15 20]*1e9
        for iter=0:7
            iter
            for alphaalpha=1
                for betabeta=0:1
                    for fault=0
                        LinearElasticity
                        
                        % % %                     B_disl=0.05;
                        % % %                     dip=1.768191886644777;
                        % % %                     trasl=4300;
                        % % %                     s1=200/cos(dip-pi/2);
                        % % %                     s2=1500/cos(dip-pi/2);
                        % % %                     [u1,v1]=InclinedDip3(Xusup-trasl,0,s1,s2,dip,ni,B_disl*2,miu);
                        % % %                     dip=1.793696559123272;
                        % % %                     trasl=7940;
                        % % %                     s1=0/cos(dip-pi/2);
                        % % %                     s2=1500/cos(dip-pi/2);
                        % % %                     [u2,v2]=InclinedDip3(Xusup-trasl,0,s1,s2,dip,ni,B_disl*2,miu);
                        % % %                     u_fault=u1+u2;
                        % % %                     v_fault=v1+v2;
                        % % %                     usup=usup+fault*u_fault';
                        % % %                     vsup=vsup+fault*v_fault';
                        
                        if h==1
                            Utime_onlyT(:,iter+1)=usup;
                            Vtime_onlyT(:,iter+1)=vsup;
                        elseif h==2
                            Utime(:,iter+1)=usup;
                            Vtime(:,iter+1)=vsup;
                            solTtime(:,iter+1)=solT(:);
                            solPtime(:,iter+1)=solP(:);
                        end
                        %                         Matrix{h}=M;
                        %                 legend_str{h}=['\mu = ',num2str(miu/1e9),' GPa'];
                        %                 legend_str{h}=['K = ',num2str(Kb/1e9),' GPa'];
                        %                         legend_str{h}=['fault = ',num2str(fault)];
                        h=h+1
                    end
                end
            end
            perc(:,iter+1)=Vtime_onlyT(:,iter+1)./Vtime(:,iter+1);
            legend_str{iter+1}=['time [yrs]  = ',num2str(time(iter+1)/31536000)];
            h=1;
        end
    end
end

figure(1)
subplot(1,2,1)
plot(Xusup,Utime)
xlim([0 10000])
xlabel('r[m]')
title('U[m]')
subplot(1,2,2)
plot(Xusup,Vtime)
xlim([0 10000])
title('V[m]')
xlabel('r[m]')
legend(num2str([1:5]'))
