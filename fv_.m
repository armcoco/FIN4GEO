function sFv = fv_(x,y)
%FV_
%    SFV = FV_(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    20-Apr-2021 14:54:11

t2 = 1.0./x;
sFv = y.*6.666666666666667e-7-t2.*(x.*5.511111111111111e2+y.*6.888888888888889e2)-t2.*(x.*2.222222222222222e-7+y.*4.444444444444444e-7).*1.24e9+6.872355555555556e5;
