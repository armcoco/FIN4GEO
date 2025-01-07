function createMatlabFunctions(n)

global AS
syms  x y xi eta;

[~,~,~,~,sG,snu,ssig,sL_SP,ssolu,ssolv,sphi1,sphi2,sxt,syt,ssolp,ssolT,salpha,sbeta,ssolgp,ssolsp]=data_symb_n(n);

matlabFunction(ssolu,'file','solu_.m','vars',[x,y]);
matlabFunction(ssolv,'file','solv_.m','vars',[x,y]);
matlabFunction(sG,'file','G_.m','vars',[x,y]);
matlabFunction(snu,'file','nu_.m','vars',[x,y]);
matlabFunction(ssig,'file','sig_.m','vars',[x,y]);
matlabFunction(sL_SP,'file','L_SP_.m','vars',[x,y]);
matlabFunction(sphi1,'file','phi1_.m','vars',[x,y]);
matlabFunction(sphi2,'file','phi2_.m','vars',[x,y]);
matlabFunction(ssolp,'file','solp_.m','vars',[x,y]);
matlabFunction(ssolT,'file','solT_.m','vars',[x,y]);
matlabFunction(salpha,'file','alpha_.m','vars',[x,y]);
matlabFunction(sbeta,'file','beta_.m','vars',[x,y]);
matlabFunction(ssolgp,'file','solgp_.m','vars',[x,y]);
matlabFunction(ssolsp,'file','solsp_.m','vars',[x,y]);

sphi1x =diff(sphi1,x,1);          sphi1y =diff(sphi1,y,1);
smodulo1 = sqrt(sphi1x^2+sphi1y^2);
snx1=sphi1x/smodulo1;               sny1=sphi1y/smodulo1;
sphi2x =diff(sphi2,x,1);          sphi2y =diff(sphi2,y,1);
smodulo2 = sqrt(sphi2x^2+sphi2y^2);
snx2=sphi2x/smodulo2;               sny2=sphi2y/smodulo2;

ssolux=diff(ssolu,x,1);           ssoluy=diff(ssolu,y,1);
ssolvx=diff(ssolv,x,1);           ssolvy=diff(ssolv,y,1);
ssolgpx=diff(ssolgp,x,1);           ssolgpy=diff(ssolgp,y,1);
ssolspx=diff(ssolsp,x,1);           ssolspy=diff(ssolsp,y,1);
matlabFunction(ssolux,'file','solux_.m','vars',[x,y]);
matlabFunction(ssoluy,'file','soluy_.m','vars',[x,y]);
matlabFunction(ssolvx,'file','solvx_.m','vars',[x,y]);
matlabFunction(ssolvy,'file','solvy_.m','vars',[x,y]);
matlabFunction(ssolgpx,'file','solgpx_.m','vars',[x,y]);
matlabFunction(ssolgpy,'file','solgpy_.m','vars',[x,y]);
matlabFunction(ssolspx,'file','solspx_.m','vars',[x,y]);
matlabFunction(ssolspy,'file','solspy_.m','vars',[x,y]);

slap_u=diff(ssolux,x,1)+diff(ssoluy,y,1);
matlabFunction(slap_u,'file','lap_u_.m','vars',[x,y]);
slap_v=diff(ssolvx,x,1)+diff(ssolvy,y,1);
matlabFunction(slap_v,'file','lap_v_.m','vars',[x,y]);

sdiv=diff(ssolu,x,1)+diff(ssolv,y,1);

sdivGgrad_u=diff(sG*ssolux,x,1)+diff(sG*ssoluy,y,1);
matlabFunction(sdivGgrad_u,'file','divGgrad_u_.m','vars',[x,y]);
sdivGgrad_v=diff(sG*ssolvx,x,1)+diff(sG*ssolvy,y,1);
matlabFunction(sdivGgrad_v,'file','divGgrad_v_.m','vars',[x,y]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sFu=sG*slap_u+1/(1-2*snu)*sG*diff(sdiv,x,1);
% matlabFunction(sFu,'file','fu_.m','vars',[x,y]);
% sFv=sG*slap_v+1/(1-2*snu)*sG*diff(sdiv,y,1);
% matlabFunction(sFv,'file','fv_.m','vars',[x,y]);
slambda=2*snu*sG/(1-2*snu);
sFu=-(diff((sG+slambda)*ssolux,x,1)+sdivGgrad_u+diff(slambda*ssolvy,x,1)+diff(sG*ssolvx,y,1));
sFv=-(diff((sG+slambda)*ssolvy,y,1)+sdivGgrad_v+diff(slambda*ssolux,y,1)+diff(sG*ssoluy,x,1));

ssigmanu1=2*sG/(1-2*snu)*((1-snu)*ssolux+snu*ssolvy)*snx1+sG*(ssoluy+ssolvx)*sny1;
ssigmanv1=sG*(ssoluy+ssolvx)*snx1+2*sG/(1-2*snu)*( snu*ssolux+(1-snu)*ssolvy)*sny1;
ssigmanu2=2*sG/(1-2*snu)*((1-snu)*ssolux+snu*ssolvy)*snx2+sG*(ssoluy+ssolvx)*sny2;
ssigmanv2=sG*(ssoluy+ssolvx)*snx2+2*sG/(1-2*snu)*( snu*ssolux+(1-snu)*ssolvy)*sny2;
% ssigmanu1=ssolu;
% ssigmanv1=ssolv;
% ssigmanu1=ssolux*snx1+ssoluy*sny1;
% ssigmanv1=ssolvx*snx1+ssolvy*sny1;
% ssigmanu2=ssolux*snx2+ssoluy*sny2;
% ssigmanv2=ssolvx*snx2+ssolvy*sny2;

if AS
    ssigmar=(slambda+2*sG)*ssolux+slambda*(ssolu/x+ssolvy)-(3*slambda+2*sG)*salpha/3*ssolT-sbeta*ssolp;
    ssigmatheta=(slambda+2*sG)*ssolu/x+slambda*(ssolux+ssolvy)-(3*slambda+2*sG)*salpha/3*ssolT-sbeta*ssolp;
    ssigmaz=(slambda+2*sG)*ssolvy+slambda*(ssolux+ssolu/x)-(3*slambda+2*sG)*salpha/3*ssolT-sbeta*ssolp;
    staurz=sG*(ssoluy+ssolvx);
    sFu=-(diff(ssigmar,x,1)+(ssigmar-ssigmatheta)/x+diff(staurz,y,1));
    sFv=-(diff(ssigmaz,y,1)+diff(staurz,x,1)+staurz/x);
    ssigmanu1=ssigmar*snx1+staurz*sny1;
    ssigmanv1=staurz*snx1+ssigmaz*sny1;
    ssigmanu2=ssigmar*snx2+staurz*sny2;
    ssigmanv2=staurz*snx2+ssigmaz*sny2;
    %%% gravitational potential
    sFgp=-(1/x*diff(x*diff(ssolgp,x,1),x,1)+diff(ssolgp,y,2));
    %%% self potential
    sJx = ssig*ssolspx+sL_SP*diff(ssolp,x,1);
    sJy = ssig*ssolspy+sL_SP*diff(ssolp,y,1);
    sFsp=-(1/x*diff(x*sJx,x,1)+diff(sJy,y,1));
end

matlabFunction(sFu,'file','fu_.m','vars',[x,y]);
matlabFunction(sFv,'file','fv_.m','vars',[x,y]);
matlabFunction(ssigmanu1,'file','sigmanu1_.m','vars',[x,y]);
matlabFunction(ssigmanv1,'file','sigmanv1_.m','vars',[x,y]);
matlabFunction(ssigmanu2,'file','sigmanu2_.m','vars',[x,y]);
matlabFunction(ssigmanv2,'file','sigmanv2_.m','vars',[x,y]);
% gravitational potential 
matlabFunction(sFgp,'file','fgp_.m','vars',[x,y]);
matlabFunction(sFsp,'file','fsp_.m','vars',[x,y]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabFunction(sxt,'file','xt_.m','vars',[xi]);
matlabFunction(syt,'file','yt_.m','vars',[eta]);
matlabFunction(diff(sxt,xi,1),'file','xt1_.m','vars',[xi]);
matlabFunction(diff(syt,eta,1),'file','yt1_.m','vars',[eta]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabFunction(snx1,sny1,'file','normal_phi1_.m','vars',[x,y],'outputs',{'nx','ny'});
matlabFunction(snx2,sny2,'file','normal_phi2_.m','vars',[x,y],'outputs',{'nx','ny'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

