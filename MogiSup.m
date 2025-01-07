function [Ux1,Uz1]=MogiSup(x,z,a,f,P0,miu,ni)
Xc=0;
Yc=0;

y=zeros(size(x));
z=zeros(size(x));
R=sqrt((x-Xc).^2+(y-Yc).^2+(z-f).^2);
Ux=P0*a^3*(x-Xc)./(R.^3);
Uy=P0*a^3*(y-Yc)./(R.^3);
Uz=P0*a^3*(z-f)./(R.^3);

Ux1=Ux*(1-ni)/(miu);
Uy1=Uy*(1-ni)/(miu);
Uz1=Uz*(1-ni)/(miu);