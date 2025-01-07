function [u,v] = GaussSeidel(u,v,dx,dy,Inside,Ghost,Mask,F,G,NU,BC,Xt1inv,Yt1inv)

% % % Ordino le celle fantasma in base alla distanza dall'interfaccia
% dists = Phi(ghost);
% [dists, I] = sort(dists);
% ghost = ghost(I);

[nn,mm]=size(u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overRel=5;
% delta=5*dx;
dtN=1*dx/4;
dtD=1*1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auLR=-2*G.*(1-NU)./(1-2*NU).*Xt1inv;
avBT=-2*G.*(1-NU)./(1-2*NU).*Yt1inv;
auBT=-G.*Yt1inv;
avLR=-G.*Xt1inv;
b_coeff1=-2*G.*NU./(1-2*NU);
c_coeff1=-G;
b_coeff=b_coeff1.*Xt1inv;
c_coeff=c_coeff1.*Yt1inv;
bc_coeff=b_coeff1+c_coeff1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter=1:overRel
    for k=Inside.B';
        gamma1=Xt1inv(k);
        gamma2=Yt1inv(k);
        auR=0.5*gamma1*(auLR(k+nn)+auLR(k));
        auL=0.5*gamma1*(auLR(k-nn)+auLR(k));
        auT=0.5*gamma2*(auBT(k+1)+auBT(k));
        auB=0.5*gamma2*(auBT(k-1)+auBT(k));
        avR=0.5*gamma1*(avLR(k+nn)+avLR(k));
        avL=0.5*gamma1*(avLR(k-nn)+avLR(k));
        avT=0.5*gamma2*(avBT(k+1)+avBT(k));
        avB=0.5*gamma2*(avBT(k-1)+avBT(k));
        bc=gamma1*gamma2*bc_coeff(k);
        bx=gamma2*Xt1inv(k)*(b_coeff1(k+nn)-b_coeff1(k-nn))/2;
        by=gamma1*Yt1inv(k)*(b_coeff1(k+1)-b_coeff1(k-1))/2;
        cx=gamma2*Xt1inv(k)*(c_coeff1(k+nn)-c_coeff1(k-nn))/2;
        cy=gamma1*Yt1inv(k)*(c_coeff1(k+1)-c_coeff1(k-1))/2;
        Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
        stencil_uu=[k k-nn k+nn k-1 k+1];
        stencil_uv=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
        stencil_vu=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
        stencil_vv=[k k-nn k+nn k-1 k+1];
        coeffs_uu=[auR+auL+auT+auB -auL -auR -auB -auT]/dx^2;
        coeffs_vv=[avR+avL+avT+avB -avL -avR -avB -avT]/dx^2;
        if all(Mask(Neigh)) %|| 1==1
            coeffs_uv=[0 cy/2 -cy/2 bx/2 -bx/2 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
            coeffs_vu=[0 by/2 -by/2 cx/2 -cx/2 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
        elseif Mask(k+nn-1) && Mask(k-nn+1)
            coeffs_uv=[bc (cy-bc)/2 -(cy+bc)/2 (bx-bc)/2 -(bx+bc)/2 0 bc/2 bc/2 0]/dx^2;
            coeffs_vu=[bc (by-bc)/2 -(by+bc)/2 (cx-bc)/2 -(cx+bc)/2 0 bc/2 bc/2 0]/dx^2;
        elseif Mask(k+nn+1) && Mask(k-nn-1)
            coeffs_uv=[-bc (cy+bc)/2 -(cy-bc)/2 (bx+bc)/2 -(bx-bc)/2 -bc/2 0 0 -bc/2]/dx^2;
            coeffs_vu=[-bc (by+bc)/2 -(by-bc)/2 (cx+bc)/2 -(cx-bc)/2 -bc/2 0 0 -bc/2]/dx^2;
        else
            disp('Lo stencil per la derivata mista include dei ghost non accettabili!')
            k
            return;
        end
        uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
        vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
        DTUinv=[coeffs_uu(1) coeffs_uv(1); coeffs_vu(1) coeffs_vv(1)];
        DTU=inv(DTUinv);
        uvk=[u(k) v(k)]'+DTU*[F.u(k)-uu F.v(k)-vv]';
        u(k)=uvk(1);
        v(k)=uvk(2);
        %     dtu=1/(coeffs_uu(1));
        %     dtv=1/(coeffs_vv(1));
        %     u(k)=u(k)+dtu*(F.u(k)-uu);
        %     v(k)=v(k)+dtv*(F.v(k)-vv);
    end
    
%     break;
    
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
        stencil_uu=current.stencil;
        stencil_uv=current.stencil_c;
        stencil_vu=current.stencil_c;
        stencil_vv=current.stencil;
        coeffs_uu=2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:);
        coeffs_uv=2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:);
        coeffs_vu=2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:);
        coeffs_vv=2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:);
        
%         coeffs_uu=current.coeffsD(:);
%         coeffs_uv=zeros(9,1);
%         coeffs_vu=zeros(9,1);
%         coeffs_vv=current.coeffsD(:);
        
        uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
        vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
        dtu=dtN/g/max(current.xt1inv,current.yt1inv);
        dtv=dtN/g/max(current.xt1inv,current.yt1inv);
        
%         dtu=dtD;
%         dtv=dtD;
        
        u(k)=u(k)+dtu*(F.u(k)-uu);
        v(k)=v(k)+dtv*(F.v(k)-vv);
    end
    
    for current=Ghost.Phi2
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
        stencil_uu=current.stencil;
        stencil_uv=current.stencil_c;
        stencil_vu=current.stencil_c;
        stencil_vv=current.stencil;
        coeffs_uu=2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:);
        coeffs_uv=2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:);
        coeffs_vu=2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:);
        coeffs_vv=2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:);
        uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
        vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
        dtu=-dtN/g/max(current.xt1inv,current.yt1inv);
        dtv=-dtN/g/max(current.xt1inv,current.yt1inv);
        u(k)=u(k)+dtu*(F.u(k)-uu);
        v(k)=v(k)+dtv*(F.v(k)-vv);
    end
    
    for current=Ghost.Bdy
        k = current.index;
        type=current.type;
        uu=u(current.stencil(:))'*current.coeffsD(:);
        vv=v(current.stencil(:))'*current.coeffsD(:);
        dtu=dtD;
        dtv=dtD;
        u(k)=u(k)+dtu*(F.u(k)-uu);
        v(k)=v(k)+dtv*(F.v(k)-vv);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=Inside.all';
    gamma1=Xt1inv(k);
    gamma2=Yt1inv(k);
    auR=0.5*gamma1*(auLR(k+nn)+auLR(k));
    auL=0.5*gamma1*(auLR(k-nn)+auLR(k));
    auT=0.5*gamma2*(auBT(k+1)+auBT(k));
    auB=0.5*gamma2*(auBT(k-1)+auBT(k));
    avR=0.5*gamma1*(avLR(k+nn)+avLR(k));
    avL=0.5*gamma1*(avLR(k-nn)+avLR(k));
    avT=0.5*gamma2*(avBT(k+1)+avBT(k));
    avB=0.5*gamma2*(avBT(k-1)+avBT(k));
    bc=gamma1*gamma2*bc_coeff(k);
    bx=gamma2*Xt1inv(k)*(b_coeff1(k+nn)-b_coeff1(k-nn))/2;
    by=gamma1*Yt1inv(k)*(b_coeff1(k+1)-b_coeff1(k-1))/2;
    cx=gamma2*Xt1inv(k)*(c_coeff1(k+nn)-c_coeff1(k-nn))/2;
    cy=gamma1*Yt1inv(k)*(c_coeff1(k+1)-c_coeff1(k-1))/2;
    Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    stencil_uu=[k k-nn k+nn k-1 k+1];
    stencil_uv=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    stencil_vu=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
    stencil_vv=[k k-nn k+nn k-1 k+1];
    coeffs_uu=[auR+auL+auT+auB -auL -auR -auB -auT]/dx^2;
    coeffs_vv=[avR+avL+avT+avB -avL -avR -avB -avT]/dx^2;
    if all(Mask(Neigh)) %|| 1==1
        coeffs_uv=[0 cy/2 -cy/2 bx/2 -bx/2 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
        coeffs_vu=[0 by/2 -by/2 cx/2 -cx/2 -bc/4 bc/4 bc/4 -bc/4]/dx^2;
    elseif Mask(k+nn-1) && Mask(k-nn+1)
        coeffs_uv=[bc (cy-bc)/2 -(cy+bc)/2 (bx-bc)/2 -(bx+bc)/2 0 bc/2 bc/2 0]/dx^2;
        coeffs_vu=[bc (by-bc)/2 -(by+bc)/2 (cx-bc)/2 -(cx+bc)/2 0 bc/2 bc/2 0]/dx^2;
    elseif Mask(k+nn+1) && Mask(k-nn-1)
        coeffs_uv=[-bc (cy+bc)/2 -(cy-bc)/2 (bx+bc)/2 -(bx-bc)/2 -bc/2 0 0 -bc/2]/dx^2;
        coeffs_vu=[-bc (by+bc)/2 -(by-bc)/2 (cx+bc)/2 -(cx-bc)/2 -bc/2 0 0 -bc/2]/dx^2;
    else
        disp('Lo stencil per la derivata mista include dei ghost non accettabili!')
        k
        return;
    end
    uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
    vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
    DTUinv=[coeffs_uu(1) coeffs_uv(1); coeffs_vu(1) coeffs_vv(1)];
    DTU=inv(DTUinv);
    uvk=[u(k) v(k)]'+DTU*[F.u(k)-uu F.v(k)-vv]';
    u(k)=uvk(1);
    v(k)=uvk(2);
    %     dtu=1/(coeffs_uu(1));
    %     dtv=1/(coeffs_vv(1));
    %     u(k)=u(k)+dtu*(F.u(k)-uu);
    %     v(k)=v(k)+dtv*(F.v(k)-vv);
end

% return;

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
    stencil_uu=current.stencil;
    stencil_uv=current.stencil_c;
    stencil_vu=current.stencil_c;
    stencil_vv=current.stencil;
    coeffs_uu=2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:);
    coeffs_uv=2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:);
    coeffs_vu=2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:);
    coeffs_vv=2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:);
    
%     coeffs_uu=current.coeffsD(:);
%     coeffs_uv=zeros(9,1);
%     coeffs_vu=zeros(9,1);
%     coeffs_vv=current.coeffsD(:);
    
    uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
    vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
    dtu=dtN/g/max(current.xt1inv,current.yt1inv);
    dtv=dtN/g/max(current.xt1inv,current.yt1inv);
    
%     dtu=dtD;
%     dtv=dtD;
    
    u(k)=u(k)+dtu*(F.u(k)-uu);
    v(k)=v(k)+dtv*(F.v(k)-vv);
end

for current=Ghost.Phi2
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
    stencil_uu=current.stencil;
    stencil_uv=current.stencil_c;
    stencil_vu=current.stencil_c;
    stencil_vv=current.stencil;
    coeffs_uu=2*g*(1-nu)/(1-2*nu)*nx*coeffs_dx(:)+g*ny*coeffs_dy(:);
    coeffs_uv=2*g*nu/(1-2*nu)*nx*coeffs_dy_c(:)+g*ny*coeffs_dx_c(:);
    coeffs_vu=2*g*nu/(1-2*nu)*ny*coeffs_dx_c(:)+g*nx*coeffs_dy_c(:);
    coeffs_vv=2*g*(1-nu)/(1-2*nu)*ny*coeffs_dy(:)+g*nx*coeffs_dx(:);
    uu=u(stencil_uu(:))'*coeffs_uu(:)+v(stencil_uv(:))'*coeffs_uv(:);
    vv=u(stencil_vu(:))'*coeffs_vu(:)+v(stencil_vv(:))'*coeffs_vv(:);
    dtu=-dtN/g/max(current.xt1inv,current.yt1inv);
    dtv=-dtN/g/max(current.xt1inv,current.yt1inv);
    u(k)=u(k)+dtu*(F.u(k)-uu);
    v(k)=v(k)+dtv*(F.v(k)-vv);
end

for current=Ghost.Bdy
    k = current.index;
    type=current.type;
    uu=u(current.stencil(:))'*current.coeffsD(:);
    vv=v(current.stencil(:))'*current.coeffsD(:);
    dtu=dtD;
    dtv=dtD;
    u(k)=u(k)+dtu*(F.u(k)-uu);
    v(k)=v(k)+dtv*(F.v(k)-vv);
end
