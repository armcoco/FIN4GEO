function snu = nu__(x,y)

snu = zeros(size(x));
%HTS
k=find(y(:)>=1500 & x(:)<50);
snu(k)=0.3;
%TZ
k=find(y(:)>=1500 & x(:)>50 & x(:)<150);
snu(k)=0.23;
%Edifice
k=find(y(:)>=1500 & x(:)>150);
snu(k)=0.17;
%Edifice bottom HTS and TZ
k=find(y(:)<1500 & y(:)>0);
snu(k)=0.17;
%Crust
k=find(y(:)<=0 & y(:)>-30000);
snu(k)=0.25;
