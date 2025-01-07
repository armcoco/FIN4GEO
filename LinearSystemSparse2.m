function [u,v,M,RHS,U]=LinearSystemSparse2(dx,dy,Inside,Ghost,Mask,F,G,NU,BC)

% w=1; %4/5; %peso dell'iterazione di Jacobi
[nn,mm]=size(F.u);

% M=speye(2*nn*mm);
RHS=nan*zeros(2*nn*mm,1);
nnmm=nn*mm;

Outside=setdiff((1:nnmm),[Inside(:);[Ghost.Phi1.index]';[Ghost.Phi2.index]';[Ghost.Bdy.index]']);
NNZ=2*numel(Outside)+14*2*numel(Inside)+18*2*numel(Ghost.Phi1)+18*2*numel(Ghost.Phi2)+3*2*numel(Ghost.Bdy);
ROWS=zeros(1,NNZ);
COLS=zeros(1,NNZ);
VALS=zeros(1,NNZ);
ROWS(1:2*numel(Outside))=[Outside,nnmm+Outside];
COLS(1:2*numel(Outside))=[Outside,nnmm+Outside];
VALS(1:2*numel(Outside))=1;
cont=2*numel(Outside);
for k=Inside';
    a11=-G(k)*2*(1-NU(k))/(1-2*NU(k));
    a22=-G(k);
    a12=-G(k)/(2-4*NU(k));
    Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    if all(Mask(Neigh)) %|| 1==1
        ROWS(cont+(1:14))=k;
        COLS(cont+(1:14))=[ [k k-nn k+nn k-1 k+1], nnmm+[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a11 -a11 -a22 -a22 0 0 0 0 0 -a12/2 a12/2 a12/2 -a12/2]/dx^2;
%         VALS(cont+(1:9))=[-4 1 1 1 1 0 0 0 0]/dx^2;
        cont=cont+14;
        RHS(k)=F.u(k);
        ROWS(cont+(1:14))=nnmm+k;
        COLS(cont+(1:14))=[nnmm+[k k-nn k+nn k-1 k+1], [k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a22 -a22 -a11 -a11 0 0 0 0 0 -a12/2 a12/2 a12/2 -a12/2]/dx^2;
%         VALS(cont+(1:9))=[-4 1 1 1 1 0 0 0 0]/dx^2;
        cont=cont+14;
        RHS(nnmm+k)=F.v(k);
    elseif Mask(k+nn-1) && Mask(k-nn+1)
        ROWS(cont+(1:14))=k;
        COLS(cont+(1:14))=[ [k k-nn k+nn k-1 k+1], nnmm+[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a11 -a11 -a22 -a22 2*a12 -a12 -a12 -a12 -a12 0 a12 a12 0]/dx^2;
%         VALS(cont+(1:14))=[2*(a11+a22+a12) -a11-a12 -a11-a12 -a22-a12 -a22-a12 0 a12 a12 0]/dx^2;
        cont=cont+14;
        RHS(k)=F.u(k);
        ROWS(cont+(1:14))=nnmm+k;
        COLS(cont+(1:14))=[nnmm+[k k-nn k+nn k-1 k+1], [k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a22 -a22 -a11 -a11 2*a12 -a12 -a12 -a12 -a12 0 a12 a12 0]/dx^2;
%         VALS(cont+(1:14))=[2*(a11+a22+a12) -a11-a12 -a11-a12 -a22-a12 -a22-a12 0 a12 a12 0]/dx^2;
        cont=cont+14;
        RHS(nnmm+k)=F.v(k);
    elseif Mask(k+nn+1) && Mask(k-nn-1)
        ROWS(cont+(1:14))=k;
        COLS(cont+(1:14))=[ [k k-nn k+nn k-1 k+1], nnmm+[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a11 -a11 -a22 -a22 -2*a12 a12 a12 a12 a12 -a12 0 0 -a12]/dx^2;
%         VALS(cont+(1:14))=[2*(a11+a22-a12) -a11+a12 -a11+a12 -a22+a12 -a22+a12 -a12 0 0 -a12]/dx^2;
        cont=cont+14;
        RHS(k)=F.u(k);
        ROWS(cont+(1:14))=nnmm+k;
        COLS(cont+(1:14))=[nnmm+[k k-nn k+nn k-1 k+1], [k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1] ];
        VALS(cont+(1:14))=[2*(a11+a22) -a22 -a22 -a11 -a11 -2*a12 a12 a12 a12 a12 -a12 0 0 -a12]/dx^2;
%         VALS(cont+(1:14))=[2*(a11+a22-a12) -a11+a12 -a11+a12 -a22+a12 -a22+a12 -a12 0 0 -a12]/dx^2;
        cont=cont+14;
        RHS(nnmm+k)=F.v(k);
    else
        disp('Lo stencil per la derivata mista include dei ghost non accettabili!')
        k
        return;
    end
end

for current=Ghost.Phi1
    k = current.index;
    nx=current.nx;
    ny=current.ny;
    g=current.G;
    nu=current.nu;
    ROWS(cont+(1:18))=k;
    COLS(cont+(1:18))=[current.stencil(:);nnmm+current.stencil_c(:)];
    VALS(cont+(1:18))=[2*g*(1-nu)/(1-2*nu)*nx*current.coeffs_dx(:)+g*ny*current.coeffs_dy(:); 2*g*nu/(1-2*nu)*nx*current.coeffs_dy_c(:)+g*ny*current.coeffs_dx_c(:)];
%     VALS(cont+(1:18))=[current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:); zeros(9,1)];
%     VALS(cont+(1:18))=[current.coeffsD(:); zeros(9,1)];
    cont=cont+18;    
    RHS(k)=F.u(k);
    ROWS(cont+(1:18))=nnmm+k;
    COLS(cont+(1:18))=[current.stencil_c(:);nnmm+current.stencil(:)];
    VALS(cont+(1:18))=[2*g*nu/(1-2*nu)*ny*current.coeffs_dx_c(:)+g*nx*current.coeffs_dy_c(:); 2*g*(1-nu)/(1-2*nu)*ny*current.coeffs_dy(:)+g*nx*current.coeffs_dx(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.coeffsD(:)];
    cont=cont+18;
    RHS(nnmm+k)=F.v(k);
end

for current=Ghost.Phi2
    k = current.index;
    nx=current.nx;
    ny=current.ny;
    g=current.G;
    nu=current.nu;    
    ROWS(cont+(1:18))=k;
    COLS(cont+(1:18))=[current.stencil(:);nnmm+current.stencil_c(:)];
    VALS(cont+(1:18))=[2*g*(1-nu)/(1-2*nu)*nx*current.coeffs_dx(:)+g*ny*current.coeffs_dy(:); 2*g*nu/(1-2*nu)*nx*current.coeffs_dy_c(:)+g*ny*current.coeffs_dx_c(:)];
%     VALS(cont+(1:18))=[current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:); zeros(9,1)];
%     VALS(cont+(1:18))=[current.coeffsD(:); zeros(9,1)];    
    cont=cont+18;    
    RHS(k)=F.u(k);
    ROWS(cont+(1:18))=nnmm+k;
    COLS(cont+(1:18))=[current.stencil_c(:);nnmm+current.stencil(:)];
    VALS(cont+(1:18))=[2*g*nu/(1-2*nu)*ny*current.coeffs_dx_c(:)+g*nx*current.coeffs_dy_c(:); 2*g*(1-nu)/(1-2*nu)*ny*current.coeffs_dy(:)+g*nx*current.coeffs_dx(:)];    
%     VALS(cont+(1:18))=[zeros(9,1);current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:)];
%     VALS(cont+(1:18))=[zeros(9,1);current.coeffsD(:)];    
    cont=cont+18;
    RHS(nnmm+k)=F.v(k);
end

for current=Ghost.Bdy
    k = current.index;
    type=current.type;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[current.stencil(:)];
    VALS(cont+(1:3))=(type==0 || type==1 || BC)*current.coeffsD(:)+((type==2 || type==3) && ~BC)*current.coeffsN(:);
    cont=cont+3;
    ROWS(cont+(1:3))=nnmm+k;
    COLS(cont+(1:3))=[nnmm+current.stencil(:)];
    VALS(cont+(1:3))=(type==2 || type==3 || BC)*current.coeffsD(:)+((type==0 || type==1) && ~BC)*current.coeffsN(:);
    cont=cont+3;    
    RHS(k)=F.u(k);
    RHS(nnmm+k)=F.v(k);
end

M = sparse(ROWS,COLS,VALS);

U=M\RHS;

I=1:nnmm;
u=reshape(U(I),nn,mm);
v=reshape(U(nnmm+I),nn,mm);
