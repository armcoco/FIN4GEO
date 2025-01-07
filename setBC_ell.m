function [F] = setBC_ell(F,Xu,Yu)

global AS

[nn,mm]=size(F.gp);
for j=2:mm-1
    i=1;
    k = i+(j-1)*nn;
    xC = Xu(k);
    yC = Yu(k);
    F.gp(k)=solgp_(xC,yC);
    F.sp(k)=solsp_(xC,yC);
    i=nn;
    k = i+(j-1)*nn;
    xC = Xu(k);
    yC = Yu(k);
    F.gp(k)=solgp_(xC,yC);
    F.sp(k)=solsp_(xC,yC);
end

for i=2:nn-1
    j=1;
    k = i+(j-1)*nn;
    xC = Xu(k);
    yC = Yu(k);
    F.gp(k)=AS*(-solgpx_(xC,yC))+(1-AS)*solgp_(xC,yC);
    F.sp(k)=AS*(-solspx_(xC,yC))+(1-AS)*solsp_(xC,yC);
    j=mm;
    k = i+(j-1)*nn;
    xC = Xu(k);
    yC = Yu(k);
    F.gp(k)=solgp_(xC,yC);
    F.sp(k)=solsp_(xC,yC);
end
