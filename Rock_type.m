function val=Rock_type(Xcenter,Ycenter,Zcenter)
Fault1P1=[4000 -1500];
Fault1P2=[4260 -200];
Fault2P1=[7600 -1500];
Fault2P2=[7940 0];

val=cell(length(Xcenter),1);
for i=1:length(Xcenter)
    x=Xcenter(i);
    y=Ycenter(i);
    z=Zcenter(i);
    P=[x z];
    if x<200
        val{i}='Fumar';
    elseif x<800
        if z<-1500
            val{i}='TCamI';
        elseif z<-1100
            val{i}='TChao';
        elseif z<-600
            val{i}='TYell';
        else
            val{i}='TPyro';
        end
    elseif x < Fault1P1(1)+(Fault1P2(1)-Fault1P1(1))/(Fault1P2(2)-Fault1P1(2))*(z-Fault1P1(2))
        if z<-1500
            val{i}='CamIn';
        elseif z<-1100
            val{i}='Chaot';
        elseif z<-600
            val{i}='Yello';
        else
            val{i}='Pyroc';
        end
    elseif x < Fault2P1(1)+(Fault2P2(1)-Fault2P1(1))/(Fault2P2(2)-Fault2P1(2))*(z-Fault2P1(2))
        if z<-1500
            val{i}='CamIn';
        elseif z<-1000
            val{i}='Chaot';
        elseif z<-500
            val{i}='Yello';
        else
            val{i}='Pyroc';
        end
    else
        if z<-1500
            val{i}='CamIn';
        elseif z<-900
            val{i}='Chaot';
        elseif z<-400
            val{i}='Yello';
        else
            val{i}='Pyroc';
        end
    end
    if OnDislocation(P,Fault1P1,Fault1P2) && abs(phi2_disl_(x,z,Fault1P1,Fault1P2))<50
        val{i}=['F',val{i}(1:4)];
    elseif OnDislocation(P,Fault2P1,Fault2P2) && abs(phi2_disl_(x,z,Fault2P1,Fault2P2))<50
        val{i}=['F',val{i}(1:4)];
    elseif x>=7600 && z<-1400 && 1==0
            val{i}='Bas  ';
    end
end

% val{2341}='FYell';
