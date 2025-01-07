function [F] = setBC(F,Ghost,BC)

global AS

for current=Ghost.Phi1
    k=current.index;
    xC = xt_(current.xc);
    yC = yt_(current.yc);
    F.u(k)=sigmanu1_(xC,yC);
    F.v(k)=sigmanv1_(xC,yC);
end

for current=Ghost.Phi2
    k=current.index;
    xC = xt_(current.xc);
    yC = yt_(current.yc);
    F.u(k)=sigmanu2_(xC,yC);
    F.v(k)=sigmanv2_(xC,yC);
end

for current=Ghost.Bdy
    k=current.index;
    type=logical(current.type);
    xC = xt_(current.xc);
    yC = yt_(current.yc);
    %     F.u(k)=(type==0 || type==1 || BC)*solu_(xC,yC)-(type==2 && ~BC)*soluy_(xC,yC)+(type==3 && ~BC)*soluy_(xC,yC);
    %     F.v(k)=(type==2 || type==3 || BC)*solv_(xC,yC)-(type==0 && ~BC)*solvx_(xC,yC)+(type==1 && ~BC)*solvx_(xC,yC);
    F.u(k)=(~(AS && ~type))*solu_(xC,yC)+(AS && ~type)*(-solux_(xC,yC));
    F.v(k)=(~(AS && ~type))*solv_(xC,yC)+(AS && ~type)*(-solvx_(xC,yC));
end

