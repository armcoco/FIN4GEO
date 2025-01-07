function sphi1 = phi1__(x,y)

sphi = zeros(size(x));

% flat topography
sphi1 = y;

% exponential decay
sphi1 = y - 3e3*exp(-(x/1e3).^2);

% Ruapehu
[ai]=find(x(:)<500);
sphi1(ai) = y(ai) - (2440 + (2650-2440)/500*x(ai));
[ai]=find(x(:)>=500 & x(:)<1200);
sphi1(ai) = y(ai) - 2650;
[ai]=find(x(:)>=1200 & x(:)<5500);
sphi1(ai) = y(ai) - (2650 + (1500-2650)/(5500-1200)*(x(ai)-1200));
[ai]=find(x(:)>=5500);
sphi1(ai) = y(ai) - 1500;
