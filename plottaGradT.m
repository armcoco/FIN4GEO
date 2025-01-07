xxx=linspace(0,10000,100);
yyy=linspace(-1500,0,20);
[XXX,YYY]=meshgrid(xxx,yyy);
Fu=interp2(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),F.u(2:end-1,2:end-1),XXX,YYY);
Fv=interp2(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),F.v(2:end-1,2:end-1),XXX,YYY);
quiver(XXX,YYY,Fu,Fv)
hold on

