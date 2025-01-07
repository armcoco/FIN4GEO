function [Rx,Ry]=RxRy_time(t)

% from Amoruso et al. (2014)
longest_semiaxis=2256;
intermediate_semiaxis=1979;
shortest_semiaxis=99;

Rx0=(longest_semiaxis+intermediate_semiaxis)/2;
Ry0=shortest_semiaxis;

% Rx0=2000;
% Ry0=100;
expansion_x=0.0;
expansion_y=0.0; %2.4587/99;

Rx=Rx0+1*expansion_x*Rx0*(1-1/(0.2*t+1));
Ry=Ry0+1*expansion_y*Ry0*(1-1/(0.2*t+1));

Rx=Rx0+1*expansion_x*Rx0*(1-exp(-0.23*t));
Ry=Ry0+1*expansion_y*Ry0*(1-exp(-0.23*t));


if t==0
    Rx=Rx0;
    Ry=Ry0;
else
    Rx=Rx0+1*expansion_x*Rx0;
    Ry=Ry0+1*expansion_y*Ry0;
end
