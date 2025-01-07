function [tank,BC,n,m,sG,snu,ssig,sL_SP,ssolu,ssolv,sphi1,sphi2,sxt,syt,ssolp,ssolT,salpha,sbeta,ssolgp,ssolsp]=data_symb_n(n)

global f miu ni P Radius m_par c_parameter c_parameterx c_parameterz m_parx ISINF AS DomInf IfPhi2 academic_test alpha beta Grav_const TRadius Rb FileTopo porosity Ks Rx Ry PPP
c_parameter=c_parameterz;
% c_parameter=1.0e4/1000; %1.5e4; %for Gottsmann 23 Jan 2018
% c_parameter=3e3; %1.5e4; %for Gottsmann 23 Jan 2018
% c_parameterx=1e4;
scale=c_parameter;

DomInf=0;
% IfPhi2=1;
ISINF=1;
m_par=1; %for Gottsmann 23 Jan 2018
% m_par=2; %for Gottsmann 23 Jan 2018
m_parx=1;

% f=0.25*scale; %5
f=3500; %0.25*scale; %5
f=3603;
P=PPP;
miu=1.24e9;
% Radius=0.15*scale; %2
Radius=500; %0.15*scale; %2
% Rx=500;
% Ry=500;
alpha=1e-5; %coefficiente di espansione termica lineare
K=5e9; %;2*miu*(1+ni)/(1-2*ni)/3; %bulk modulus K=lambda+2/3miu=E/(1-2*ni)/3=2miu(1+ni)/(1-2ni)/3
ni=(3*K-2*miu)/(2*miu+6*K);
Ks=30e9;
beta=1*(1-K/Ks);
Grav_const=6.672*1e-11;
TRadius=50;
Rb=600;
FileTopo='..\Topografia\HalfSpace.txt';
porosity=0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tank=[AS-1 1 -1 1]*(1-0.1*academic_test); %*40*scale;
BC=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=40;
m=n*(tank(2)-tank(1))/(tank(4)-tank(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~academic_test
    return
end

syms x y;

% ssol_i = 1*((x-cc)^2+(y-cc)^2);%exp(x*y);
% ssol_e = 1*((x-cc)^2+(y-cc)^2);

ssolu = x*y+x+y;
ssolv = 4*x*y+3*x-y;

% ssol_i = exp(x+y);
% ssol_e = exp(x+y); %   exp(-(x+y));
% 
% ssol_i = 1*exp(x^2+y^2);
% ssol_e = 1*exp(x^2+y^2); %   exp(-(x+y));

x0=0.0;
y0=0.0;

if ~exist('scale','var') || isempty(scale)
    scale=1;
end
ssolu = 1/scale*((x-x0)+(y-y0));
ssolv = 2/scale*((x-x0)+(y-y0));

ssolu = 2/scale*((x-x0)^2+(y-y0)^2);
ssolv = 1/scale*((x-x0)^2+(y-y0)^2);

ssolu = 1*2*((x/scale-x0)^2+(y/scale-y0)^2+(x/scale-x0)*(y/scale-y0));
ssolv = 1*1*((x/scale-x0)^2+(y/scale-y0)^2+(x/scale-x0)*(y/scale-y0));
ssolgp =3*1*((x/scale-x0)^2+(y/scale-y0)^2);
ssolsp =3*1*((x/scale-x0)^2+(y/scale-y0)^2);
% ssolsp =exp(-((x/scale-x0)^2+(y/scale-y0)^2));

% ssolT = 1*((x/scale-x0)^2+(y/scale-y0)^2+(x/scale-x0)^2*(y/scale-y0)^2);
% ssolp = 2*((x/scale-x0)^2+(y/scale-y0)^2+(x/scale-x0)^2*(y/scale-y0)^2);
ssolT = 1*((x/scale-x0)+(y/scale-y0));
ssolp = 4*(x/scale)^2 + 3*(y/scale)^2 + 2;
salpha=1; %10*(x/scale+y/scale);
sbeta=1; %20*(x/scale+y/scale);

% ssolu = 2*exp(-((x/scale)^2+(y/scale)^2));
% ssolv = 1*exp(-((x/scale)^2+(y/scale)^2));

% ssolu = 2/scale*((x-x0)^3+(y-y0)^3+(x-x0)^2*(y-y0)^2);
% ssolv = 1/scale*((x-x0)^3+(y-y0)^3+(x-x0)^2*(y-y0)^2);

% ssolu = exp(-((x/scale)^2+(y/scale)^2));
% ssolv = exp(-((x/scale)^2+(y/scale)^2));

% ssolu = exp(-(x/scale+y/scale));
% ssolv = exp(-(x/scale+y/scale));

% ssolu = 2/scale*exp(-((x/1e4)^2+(y/1e4)^2));
% ssolv = 1/scale*exp(-((x/1e4)^2+(y/1e4)^2));

% ssolu = 2/scale*((x-x0)+(y-y0));
% ssolv = 1/scale*((x-x0)+(y-y0));
% ssolu = 2;
% ssolv = 1;

F=f;
z=-y;
A=Radius^2*P/(2*miu);
u0=A*x*(1/(x^2+(z-F)^2)+1/(x^2+(z+F)^2));
v0=A*((z-F)/(x^2+(z-F)^2)+(z+F)/(x^2+(z+F)^2));
up=-A*x/(x^2+(z+F)^2)^2*(3*z^2+2*F*z-F^2-x^2);
vp=-A*1/(x^2+(z+F)^2)^2*(5*z^3+13*F*z^2+11*F^2*z+3*F^3+3*F*x^2+x^2*z);
% ssolu = u0+up;
% ssolv = -(v0+vp);

% ssolu = 2*((x-x0)^2+(y-y0)^2);
% ssolv = 1*((x-x0)^2+(y-y0)^2);

% ss=0.1;
% ssolu = sin(ss*x/scale)*sin(ss*y/scale);
% ssolv = sin(ss*x/scale)*sin(ss*y/scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sG=miu; %*(1+10*scale*x+y);
snu=0.25; %+1e-6*(y-3*x);
ssig = (x/scale)^2 + (y/scale)^2 + 1;
sL_SP = 3*(x/scale)^2 + 2*(y/scale)^2 + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x10=0.0*sqrt(2)/30;
y10=-0*0.5*scale; % sqrt(3)/40;
x20=0.0*scale;
y20=-f; %-0.5*scale;
% x10=0.05;
% y10=0.05;
% x20=0.05;
% y20=0.05;
R1=0.5*scale;
R2=Radius; %0.6*scale;
sphi1 = sqrt((x-x10)^2+(y-y10)^2)-R1;
sphi2 = sqrt((x-x20)^2+(y-y20)^2)-R2;
% sphi2 = (x-x20)^2+(y-y20)^2-R2^2;
% sphi2 = 1;

ai=0.9*scale;
bi=0.9*scale;
xi=0.0*scale;
% sphi1=y-(ai^2*bi)/((x-xi)^2+ai^2);

sphi1 = y;
% sphi1 = y - 3e3*exp(-(x/1e3).^2);

% sphi1 = y-3*scale;
% sphi2 = x^2+y^2+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms xi eta;
sxt=-log(-xi);
syt=-log(41-eta);
sxt=10*xi/sqrt(tank(2)^2-xi^2);
syt=10*eta/sqrt(tank(4)^2-eta^2);
if ~ISINF
    sxt=c_parameter*xi;
    syt=c_parameter*eta;
else
    % sxt=c_parameter*xi/(tank(2)^2-xi^2)^m_par;
    % syt=c_parameter*eta/(tank(4)^2-eta^2)^m_par;
    sxt=c_parameter*xi/(1-xi^2)^m_par;
    syt=c_parameter*eta/(1-eta^2)^m_par;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x0=0.03*sqrt(3);
% y0=0.04*sqrt(2);
% r1=0.45;
% r2=1/7; %0.2;
% sphi1 = sqrt((x-x0).^2+(y-y0).^2)-r1-((y-y0).^5+5*(x-x0).^4.*(y-y0)-10*(x-x0).^2.*(y-y0).^3)*r2./(sqrt((x-x0).^2+(y-y0).^2).^5);
% 
% x0=0.03*sqrt(3);
% y0=0.04*sqrt(2);
% r1=0.75;
% r2=1/8; %0.2;
% sphi2 = sqrt((x-x0).^2+(y-y0).^2)-r1-((y-y0).^5+5*(x-x0).^4.*(y-y0)-10*(x-x0).^2.*(y-y0).^3)*r2./(sqrt((x-x0).^2+(y-y0).^2).^5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

