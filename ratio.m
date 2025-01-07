a=4.92;
b=0.14;
c=30198;

dV_magma=a*(1-exp(-b*15))*1e9;
dM_magma=dV_magma*2350;
dM_fluid=6050-3400+c+61/60*c;
dM_fluid=dM_fluid*365.25*15*1e3;
ratio=dM_fluid/dM_magma;