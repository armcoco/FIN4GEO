function [val,RockDensity]=Por(Xcenter,Ycenter,Zcenter)
Fault1P1=[4000 -1500];
Fault1P2=[4260 -200];
Fault2P1=[7600 -1500];
Fault2P2=[7940 0];

% % % Fault1P1=[1500 -3000];
% % % Fault1P2=[1760 -200];
% % % Fault2P1=[7600 -3000];
% % % Fault2P2=[7940 0];

val=zeros(size(Xcenter));
RockDensity=zeros(size(Xcenter));
for i=1:length(Xcenter)
    x=Xcenter(i);
    y=Ycenter(i);
    z=Zcenter(i);
    P=[x z];
    if x<200
        val(i)=0.1;
        RockDensity(i)=1800;
    elseif x<800
        if z<-1500
            val(i)=0.1;
            RockDensity(i)=2000;
        elseif z<-1100
            val(i)=0.15;
            RockDensity(i)=2000;
        elseif z<-600
            val(i)=0.15;
            RockDensity(i)=1800;
        else
            val(i)=0.15;
            RockDensity(i)=1600;
        end
    elseif x < Fault1P1(1)+(Fault1P2(1)-Fault1P1(1))/(Fault1P2(2)-Fault1P1(2))*(z-Fault1P1(2))
        if z<-1500
            val(i)=0.1;
            RockDensity(i)=2000;
        elseif z<-1100
            val(i)=0.15;
            RockDensity(i)=2000;
        elseif z<-600
            val(i)=0.35;
            RockDensity(i)=1800;
        else
            val(i)=0.45;
            RockDensity(i)=1600;
        end
    elseif x < Fault2P1(1)+(Fault2P2(1)-Fault2P1(1))/(Fault2P2(2)-Fault2P1(2))*(z-Fault2P1(2))
        if z<-1500
            val(i)=0.1;
            RockDensity(i)=2000;
        elseif z<-1000
            val(i)=0.15;
            RockDensity(i)=2000;
        elseif z<-500
            val(i)=0.35;
            RockDensity(i)=1800;
        else
            val(i)=0.45;
            RockDensity(i)=1600;
        end
    else
        if z<-1500
            val(i)=0.1;
            RockDensity(i)=2000;
        elseif z<-900
            val(i)=0.15;
            RockDensity(i)=1800;
        elseif z<-400
            val(i)=0.35;
            RockDensity(i)=1600;
        else
            val(i)=0.45;
            RockDensity(i)=1600;
        end
    end
    if x>=7600 && z<-1400 && 1==0
            val(i)=0.1;
            RockDensity(i)=2000;
    end
end

% val{2341}='FYell';
