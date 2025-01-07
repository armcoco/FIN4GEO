function [u,v,M,RHS,U]=LinearSystem(dx,dy,Inside,Ghost,Mask,F,G,NU,BC)

% w=1; %4/5; %peso dell'iterazione di Jacobi
[nn,mm]=size(F.u);

M=speye(2*nn*mm);
RHS=zeros(2*nn*mm,1);

nnmm=nn*mm;

for k=Inside';
    a11=-G(k)*2*(1-NU(k))/(1-2*NU(k));
    a22=-G(k);
    a12=-G(k)/(2-4*NU(k));
    Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    if all(Mask(Neigh))
        M(k,[ [k k-nn k+nn k-1 k+1], nnmm+[k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
            [2*(a11+a22) -a11 -a11 -a22 -a22 -a12/2 a12/2 a12/2 -a12/2]/dx^2;
        RHS(k)=F.u(k);        
        M(nnmm+k,[nnmm+[k k-nn k+nn k-1 k+1], [k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
            [2*(a11+a22) -a22 -a22 -a11 -a11 -a12/2 a12/2 a12/2 -a12/2]/dx^2;
        RHS(nnmm+k)=F.v(k);
    elseif Mask(k+nn-1) && Mask(k-nn+1)
            M(k,[ [k k-nn k+nn k-1 k+1], nnmm+[k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
                [2*(a11+a22+a12) -a11-a12 -a11-a12 -a22-a12 -a22-a12 0 a12 a12 0]/dx^2;
            RHS(k)=F.u(k);
            M(nnmm+k,[nnmm+[k k-nn k+nn k-1 k+1], [k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
                [2*(a11+a22+a12) -a11-a12 -a11-a12 -a22-a12 -a22-a12 0 a12 a12 0]/dx^2;
            RHS(nnmm+k)=F.v(k);            
    elseif Mask(k+nn+1) && Mask(k-nn-1)
            M(k,[ [k k-nn k+nn k-1 k+1], nnmm+[k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
                [2*(a11+a22-a12) -a11+a12 -a11+a12 -a22+a12 -a22+a12 -a12 0 0 -a12]/dx^2;
            RHS(k)=F.u(k);
            M(nnmm+k,[nnmm+[k k-nn k+nn k-1 k+1], [k+nn+1 k+nn-1 k-nn+1 k-nn-1] ])=...
                [2*(a11+a22-a12) -a11+a12 -a11+a12 -a22+a12 -a22+a12 -a12 0 0 -a12]/dx^2;
            RHS(nnmm+k)=F.v(k);          
    else
        disp('Lo stencil per la derivata mista include dei ghost non accettabili!')
        return;
    end
end

for current=Ghost.Phi1
    k = current.index;
    M(k,[current.stencil(:)])=current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:);
    RHS(k)=F.u(k);
    M(nnmm+k,nnmm+[current.stencil(:)])=current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:);
    RHS(nnmm+k)=F.v(k);
end

for current=Ghost.Phi2
    k = current.index;
    M(k,[current.stencil(:)])=current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:);
    RHS(k)=F.u(k);
    M(nnmm+k,nnmm+[current.stencil(:)])=current.nx*current.coeffs_dx(:)+current.ny*current.coeffs_dy(:);
    RHS(nnmm+k)=F.v(k);
end

for current=Ghost.Bdy
    k = current.index;
    type=current.type;
    M(k,[current.stencil(:)])=(type==0 || type==1 || BC)*current.coeffsD(:)+((type==2 || type==3) && ~BC)*current.coeffsN(:);
    M(nnmm+k,nnmm+[current.stencil(:)])=(type==2 || type==3 || BC)*current.coeffsD(:)+((type==0 || type==1) && ~BC)*current.coeffsN(:);
    RHS(k)=F.u(k);
    RHS(nnmm+k)=F.v(k);
end

U=M\RHS;

I=1:nnmm;
u=reshape(U(I),nn,mm);
v=reshape(U(nnmm+I),nn,mm);
