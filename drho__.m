function val=drho__(x,y)

val=poros__(x,y).*(rhof__(x,y)-rhog__(x,y)).*(satur__(x,y)-1);


global f Radius;

x0=0;
y0=-f;

sphi2 = sqrt((x-x0).^2+(y-y0).^2)-Radius;

val=200*ones(size(x));
val(find(sphi2>0))=0;