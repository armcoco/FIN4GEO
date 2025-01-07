function val=phi2_disl_(x,y,P1,P2)

x1=P1(1);
y1=P1(2);
x2=P2(1);
y2=P2(2);
if x1==x2
    val=-((y1-y2)*x+x1*y2-y1*x2)/sqrt((y1-y2)^2);
else
    val=-((y1-y2)*x+(x2-x1)*y+x1*y2-y1*x2)/sqrt((y1-y2)^2+(x2-x1)^2);
%     val=((y-y2)*(x1-x2)-(x-x2)*(y1-y2));
end
