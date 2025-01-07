function [u,v,M,RHS,U,Pinv,B]=LinearSystemSparse2het(dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv)

global DomInf IsPhi2

% w=1; %4/5; %peso dell'iterazione di Jacobi
[nn,mm]=size(F.u);

% M=speye(2*nn*mm);
RHS=0*zeros(2*nn*mm,1);
Pinv_diag=zeros(2*nn*mm,1);
nnmm=nn*mm;

auLR=2*G.*(1-NU)./(1-2*NU).*Xt1inv;
avBT=2*G.*(1-NU)./(1-2*NU).*Yt1inv;
auBT=G.*Yt1inv;
avLR=G.*Xt1inv;
b_coeff1=2*G.*NU./(1-2*NU);
c_coeff1=G;
bc_coeff=b_coeff1+c_coeff1;

% auLR=-Xt1inv;
% auBT=-Yt1inv;
% avLR=-Xt1inv;
% avBT=-Yt1inv;
% b_coeff=zeros(nn,mm);
% c_coeff=zeros(nn,mm);
% bc_coeff=zeros(nn,mm);

if DomInf
    if IsPhi2 Outside=setdiff((1:nnmm),[Inside.all(:);[Ghost.Phi2.index]';[Ghost.Bdy.index]']);
    else Outside=setdiff((1:nnmm),[Inside.all(:);[Ghost.Bdy.index]']);
    end
else
    if IsPhi2 Outside=setdiff((1:nnmm),[Inside.all(:);[Ghost.Phi1.index]';[Ghost.Phi2.index]';[Ghost.Bdy.index]']);
    else Outside=setdiff((1:nnmm),[Inside.all(:);[Ghost.Phi1.index]';[Ghost.Bdy.index]']);
    end
end

% Outside=setdiff((1:nnmm),[Inside.all(:);[Ghost.Phi1.index]';[Ghost.Phi2.index]';[Ghost.Bdy.index]']);
NNZ=2*numel(Outside)+14*2*numel(Inside.all)+18*2*numel(Ghost.Phi1)+18*2*numel(Ghost.Phi2)+3*2*numel(Ghost.Bdy);
ROWS=zeros(1,NNZ);
COLS=zeros(1,NNZ);
VALS=zeros(1,NNZ);
ROWS(1:2*numel(Outside))=[Outside,nnmm+Outside];
COLS(1:2*numel(Outside))=[Outside,nnmm+Outside];
VALS(1:2*numel(Outside))=1;
cont=2*numel(Outside);
for k=Inside.all';
    gamma1=Xt1inv(k);
    gamma2=Yt1inv(k);
    gamma12=gamma1*gamma2;
    auR=0.5*gamma1*(auLR(k+nn)+auLR(k));
    auL=0.5*gamma1*(auLR(k-nn)+auLR(k));
    auT=0.5*gamma2*(auBT(k+1)+auBT(k));
    auB=0.5*gamma2*(auBT(k-1)+auBT(k));
    avR=0.5*gamma1*(avLR(k+nn)+avLR(k));
    avL=0.5*gamma1*(avLR(k-nn)+avLR(k));
    avT=0.5*gamma2*(avBT(k+1)+avBT(k));
    avB=0.5*gamma2*(avBT(k-1)+avBT(k));
    bc=gamma12*bc_coeff(k);
    bx=gamma12*(b_coeff1(k+nn)-b_coeff1(k-nn))/2;
    by=gamma12*(b_coeff1(k+1)-b_coeff1(k-1))/2;
    cx=gamma12*(c_coeff1(k+nn)-c_coeff1(k-nn))/2;
    cy=gamma12*(c_coeff1(k+1)-c_coeff1(k-1))/2;
    Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    VALS1=[auR+auL+auT+auB -auL -auR -auB -auT 0 cy/2 -cy/2 bx/2 -bx/2 0 0 0 0]/dx^2;
    VALS2=[avR+avL+avT+avB -avL -avR -avB -avT 0 by/2 -by/2 cx/2 -cx/2 0 0 0 0]/dx^2;
    if all(Mask(Neigh)) %|| 1==1
        VALS1=VALS1+[0 0 0 0 0 0 0 0 0 0 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
        VALS2=VALS2+[0 0 0 0 0 0 0 0 0 0 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
    elseif Mask(k+nn-1) && Mask(k-nn+1)
        VALS1=VALS1+[0 0 0 0 0 bc -bc/2 -bc/2 -bc/2 -bc/2 0 bc/2 bc/2 0]/dx^2;
        VALS2=VALS2+[0 0 0 0 0 bc -bc/2 -bc/2 -bc/2 -bc/2 0 bc/2 bc/2 0]/dx^2;
    elseif Mask(k+nn+1) && Mask(k-nn-1)
        VALS1=VALS1+[0 0 0 0 0 -bc bc/2 bc/2 bc/2 bc/2 -bc/2 0 0 -bc/2]/dx^2;
        VALS2=VALS2+[0 0 0 0 0 -bc bc/2 bc/2 bc/2 bc/2 -bc/2 0 0 -bc/2]/dx^2;
    else
        disp('Lo stencil per la derivata mista include dei ghost non accettabili!')
        k
        return;
    end
    ROWS(cont+(1:14))=k;
    COLS(cont+(1:14))=[ [k k-nn k+nn k-1 k+1], nnmm+[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
    VALS(cont+(1:14))=VALS1;
    cont=cont+14;
    RHS(k)=F.u(k);
    ROWS(cont+(1:14))=nnmm+k;
    COLS(cont+(1:14))=[nnmm+[k k-nn k+nn k-1 k+1], [k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
    VALS(cont+(1:14))=VALS2;
    cont=cont+14;
    RHS(nnmm+k)=F.v(k);
    Pinv_diag(k)=dx^2/(auR+auL+auT+auB);
    Pinv_diag(nnmm+k)=dx^2/(avR+avL+avT+avB);
end

for current=Ghost.Phi1
    k = current.index;
    nx=current.xt1inv*current.nx;
    ny=current.yt1inv*current.ny;
    modl=sqrt(nx^2+ny^2);
    nx=nx/modl;
    ny=ny/modl;
%     nx=current.nx;
%     ny=current.ny;
    g=current.G;
    nu=current.nu;
    coeffs_dx=current.xt1inv*current.coeffs_dx;
    coeffs_dy=current.yt1inv*current.coeffs_dy;
    coeffs_dx_c=current.xt1inv*current.coeffs_dx_c;
    coeffs_dy_c=current.yt1inv*current.coeffs_dy_c;    
    ROWS(cont+(1:18))=k;
    COLS(cont+(1:18))=[current.stencil(:);nnmm+current.stencil_c(:)];
    VALS(cont+(1:18))=[2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:); 2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:)];
%     VALS(cont+(1:18))=[current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:); zeros(9,1)];
%     VALS(cont+(1:18))=[current.coeffsD(:); zeros(9,1)];
%     VALS(cont+(1:18))=[coeffs_dy(:); zeros(9,1)];
    cont=cont+18;    
    RHS(k)=F.u(k);
    ROWS(cont+(1:18))=nnmm+k;
    COLS(cont+(1:18))=[current.stencil_c(:);nnmm+current.stencil(:)];
    VALS(cont+(1:18))=[2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:); 2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.coeffsD(:)];
    cont=cont+18;
    RHS(nnmm+k)=F.v(k);
    Pinv_diag(k)=dx/8/g/max(current.xt1inv,current.yt1inv);
    Pinv_diag(nnmm+k)=dx/8/g/max(current.xt1inv,current.yt1inv);
end

for current=Ghost.Phi2
    k = current.index;
    nx=current.xt1inv*current.nx;
    ny=current.yt1inv*current.ny;
    modl=sqrt(nx^2+ny^2);
    nx=nx/modl;
    ny=ny/modl;
    g=current.G;
    nu=current.nu;    
    coeffs_dx=current.xt1inv*current.coeffs_dx;
    coeffs_dy=current.yt1inv*current.coeffs_dy;
    coeffs_dx_c=current.xt1inv*current.coeffs_dx_c;
    coeffs_dy_c=current.yt1inv*current.coeffs_dy_c;    
    ROWS(cont+(1:18))=k;
    COLS(cont+(1:18))=[current.stencil(:);nnmm+current.stencil_c(:)];
    VALS(cont+(1:18))=[2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:); 2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:)];
%     VALS(cont+(1:18))=[current.nx*coeffs_dx(:)+current.ny*coeffs_dy(:);zeros(9,1)];
%     VALS(cont+(1:18))=[current.coeffsD(:); zeros(9,1)];
%     VALS(cont+(1:18))=[nx*coeffs_dx(:)+ny*coeffs_dy(:); zeros(9,1)];
    cont=cont+18;    
    RHS(k)=F.u(k);
    ROWS(cont+(1:18))=nnmm+k;
    COLS(cont+(1:18))=[current.stencil_c(:);nnmm+current.stencil(:)];
    VALS(cont+(1:18))=[2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:); 2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.nx*coeffs_dx(:)+current.ny*coeffs_dy(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.coeffsD(:)];    
    cont=cont+18;
    RHS(nnmm+k)=F.v(k);
    Pinv_diag(k)=-dx/8/g/max(current.xt1inv,current.yt1inv);
    Pinv_diag(nnmm+k)=-dx/8/g/max(current.xt1inv,current.yt1inv);    
end

for current=Ghost.Bdy
    k = current.index;
    type=current.type;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[current.stencil(:)];
    VALS(cont+(1:3))=current.coeffsD(:);
    cont=cont+3;
    ROWS(cont+(1:3))=nnmm+k;
    COLS(cont+(1:3))=[nnmm+current.stencil(:)];
    VALS(cont+(1:3))=current.coeffsD(:);
    cont=cont+3;    
    RHS(k)=F.u(k);
    RHS(nnmm+k)=F.v(k);
    Pinv_diag(k)=0.8;
    Pinv_diag(nnmm+k)=0.8;
end

M = sparse(ROWS,COLS,VALS);

U=M\RHS;

I=1:nnmm;
u=reshape(U(I),nn,mm);
v=reshape(U(nnmm+I),nn,mm);

Pinv=diag(Pinv_diag);
B=speye(2*nnmm)-Pinv*M;