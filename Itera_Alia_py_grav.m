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
                for betabeta=1
                    for fault=0
                        LinearElasticity
                        dGtime(:,iter+1)=dGsupT*1e8;
                        h=h+1
                    end
                end
            end
            legend_str{iter+1}=['time [yrs]  = ',num2str(time(iter+1)/31536000)];
            h=1;
        end
    end
end

figure(1)
subplot(1,2,1)
plot(traccia(:,1),dGtime)
xlim([0 3000])
xlabel('r[m]')
title('dG[miuGal]')
% legend(num2str([0.1:0.1:3]'))
legend(legend_str)

% figure(2)
subplot(1,2,2)
plot([0.1:0.1:3]',dGtime(2,:))
xlabel('time[y]')
title('dG[miuGal]')
xlim([0 3])
% legend(legend_str)
