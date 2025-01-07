function ssolT = solT__(x,y)
global TRadius Rb Radius f
ssolT = zeros(size(x));
%ssolT = 0*(x+y);
Ra=Radius;
Ta=TRadius;

r=sqrt(x.^2+(y+f).^2);

%Linear Function in a shell Radius and Rb
At=Ta/(Ra-Rb);
Bt=-Ta*Rb/(Ra-Rb);
ssolT=(At.*r+Bt).*(r<Rb);

%Parabolic Function
ssolT=Ta*1e4./r.^2;

ssolT(1,:)=0;
ssolT(end,:)=0;
ssolT(:,end)=0;
%ssolT(129,1)=0;

ssolT=Ta.*(r<Radius);
