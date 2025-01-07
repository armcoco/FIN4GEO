function [por,rho]=por_and_rho(Xcenter,Ycenter,Zcenter)

por=zeros(size(Xcenter));
rho=zeros(size(Xcenter));

for i=1:length(Xcenter)
    x=Xcenter(i);
    y=Ycenter(i);
    z=Zcenter(i);
    P=[x z];
    if x<50
        por(i)=0.20;
        rho(i)=2200;
    elseif x>50 && x<150 
        por(i)=0.10;
        rho(i)=2400;
    else
        por(i)=0.05;
        rho(i)=2500;   
    end
end
