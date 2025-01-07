function sG = G__(x,y)

sG = zeros(size(x));

%HTS
k=find(y(:)>=1500 & x(:)<50);
sG(k)=8.1301e9;
%TZ
k=find(y(:)>=1500& x(:)>50 &x(:)<150);
sG(k)=8.1301E9;
%Edifice
k=find(y(:)>=1500 & x(:)>150 );
sG(k)=1.2821e10;
%Edifice 2
k=find(y(:)<1500);
sG(k)=1.2821e10;
% Crust
% k=find(y(:)<0);
% sG(k)=0.0005*y^3+0.0162*y^2-0.2462*y+12.998;


