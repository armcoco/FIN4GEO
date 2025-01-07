function sphi2 = phi2__(x,y)

global f Rx Ry;

x0=0;
y0=-f;

sphi2 = sqrt(((x-x0)/Rx).^2+((y-y0)/Ry).^2)-1;
