function Gravity=GravSol(x,z,a,f,dRho)
Xc=0;
Yc=0;
G=6.672*1e-11;
y=zeros(size(x));
% z=zeros(size(x));
R=sqrt((x-Xc).^2+(y-Yc).^2+(z-f).^2);
Gravity=1e8*G*dRho*4/3*pi*a^3*(z-f)./(R.^3);
Gravity=1e8*G*2.4*1e6*(z-f)./(R.^3);
