clear all
close all

t = 0;
T = 1;

a = 0; b = 1;
N = 100;

dx = (b-a)/N;

x = a+dx/2:dx:b-dx/2;
x = x( : );

% Initial Condition
y = 2+t*cos(2*pi*x);

again = 1;
Id = speye ( N ) ;

e = ones(N,1);
Lap = spdiags([e -2*e e],[-1:1],N,N);
Lap(1,N) = 1;
Lap(N,1) = 1;
dt = 0.01*T/N;
% dt=dx^2;
while (t<T && again)
    
    if t+dt>T
        again = 0;
        dt = T-t;
    end
    fun = cos(2*pi*x)+t*(2+t*cos(2*pi*x))*4*pi^2.*cos(2*pi*x);
    y_star = y + dt*fun;
    ys = spdiags(y,0,N,N);
%     y=y_star+dt*ys*Lap/dx^2*y;
    y = (Id-dt*ys*Lap/dx^2)\y_star;
    
%     y=y-mean(y)+1;
    t = t+dt;
    yexa = 2+t*cos(2*pi*x);
    plot(x,y,x,yexa,'--');
    drawnow
    t
    
end