function [Fu,Fv]=GradSource(Xu,Yu,depth,radius,ampl)

Fu=Xu;
Fv=Yu+depth;
modl=sqrt(Fu.^2+Fv.^2);
Fu=Fu./modl*ampl;
Fv=Fv./modl*ampl;

Fu(sqrt(Xu.^2+(Yu+depth).^2)>radius)=0;
Fv(sqrt(Xu.^2+(Yu+depth).^2)>radius)=0;