function ssigmanv1 = sigmanv1_(x,y)
%SIGMANV1_
%    SSIGMANV1 = SIGMANV1_(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Apr-2021 14:54:12

t2 = x.^2;
t3 = y.^2;
ssigmanv1 = t2.*(-4.444444444444444e-7)-t3.*3.333333333333333e-7-x.*6.879244444444445e5-y.*6.877866666666667e5+(t2.*2.755555555555556e2+t3.*2.755555555555556e2+x.*y.*2.755555555555556e2)./x-2.0;
