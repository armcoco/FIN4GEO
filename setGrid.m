function [Inside,Ghost,Mask,MaskI]=setGrid(Phi1,Phi2,x,y,dx,dy,G,NU,Nx,Ny,Xt1inv,Yt1inv)

global FLAG c_parameter;
s=3;
[nn,mm]=size(Phi1);

delta=5*dx;
deltau=delta*c_parameter;

%%% setting the Inside structure
[X,Y]=meshgrid(x,y);
Inside.all = find(Phi1<0 & Phi2>=0 & X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end));
% Inside.B = find(Phi1<0 & Phi2>=0 & (Phi1>-deltau | Phi2<deltau) & X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end));
Inside.B = find(Phi1<0 & Phi2>=0 & X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end) & (Phi1>-deltau | Phi2<deltau | abs(X-X(1))<delta | abs(X-X(end))<delta | abs(Y-Y(1))<delta | abs(Y-Y(end))<delta ));

%%% setting the Ghost structure
Ghost=struct('Phi1',[],'Phi2',[],'Bdy',[]);
IsGhost=struct('Phi1',zeros(nn,mm),'Phi2',zeros(nn,mm),'Bdy',zeros(nn,mm));

for j=2:mm-1
    for i=2:nn-1
        k=i+(j-1)*nn;
        Neigh=[k+nn,k-nn,k+1,k-1];
%         Neigh=[k k-nn k+nn k-1 k+1 k+nn+1 k+nn-1 k-nn+1 k-nn-1];
         NeighCorner1=[k+nn+1 k-nn-1];
         NeighCorner2=[k-nn+1 k+nn-1];
        if Phi1(k)>=0
            if any(Phi1(Neigh)<0) || any(Phi1(NeighCorner1)<0 & Nx.Phi1(NeighCorner1).*Ny.Phi1(NeighCorner1)<0) || any(Phi1(NeighCorner2)<0 & Nx.Phi1(NeighCorner2).*Ny.Phi1(NeighCorner2)>=0)
                Ghost.Phi1(end+1).index=k;
                IsGhost.Phi1(k)=1;
            end
        end
        if Phi2(k)<0
            if any(Phi2(Neigh)>=0) || any(Phi2(NeighCorner1)>=0 & Nx.Phi2(NeighCorner1).*Ny.Phi2(NeighCorner1)<0) || any(Phi2(NeighCorner2)>=0 & Nx.Phi2(NeighCorner2).*Ny.Phi2(NeighCorner2)>=0)
                Ghost.Phi2(end+1).index=k;
                IsGhost.Phi2(k)=1;
            end
        end
    end
end

for j=2:mm-1
    i=1;
    k=i+(j-1)*nn;
    if Phi1(k+1)<0 && Phi2(k+1)>=0 %|| Phi1(k+1+nn)<0 || Phi1(k+1-nn)<0
        Ghost.Bdy(end+1).index=k;
        Ghost.Bdy(end).type=2;
        IsGhost.Bdy(k)=1;
    end
    i=nn;
    k=i+(j-1)*nn;
    if Phi1(k-1)<0 && Phi2(k-1)>=0 %|| Phi1(k-1+nn)<0 || Phi1(k-1-nn)<0
        Ghost.Bdy(end+1).index=k;
        Ghost.Bdy(end).type=3;
        IsGhost.Bdy(k)=1;
    end    
end

for i=2:nn-1
    j=1;
    k=i+(j-1)*nn;
    if Phi1(k+nn)<0 && Phi2(k+nn)>=0 %|| Phi1(k+nn+1)<0 || Phi1(k+nn-1)<0
        Ghost.Bdy(end+1).index=k;
        Ghost.Bdy(end).type=0;
        IsGhost.Bdy(k)=1;
    end
    j=mm;
    k=i+(j-1)*nn;
    if Phi1(k-nn)<0 && Phi2(k-nn)>=0 %|| Phi1(k-nn+1)<0 || Phi1(k-nn-1)<0
        Ghost.Bdy(end+1).index=k;
        Ghost.Bdy(end).type=1;
        IsGhost.Bdy(k)=1;
    end
end

Mask=( (Phi1<0 & Phi2>=0 & X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end)) | IsGhost.Phi1 | IsGhost.Phi2 | IsGhost.Bdy);
MaskI=( (Phi1<0 & Phi2>=0 & X>X(1) & X<X(end) & Y>Y(1) & Y<Y(end)));

for i=1:length(Ghost.Phi1)
    FLAG=1;
    k=Ghost.Phi1(i).index;
%     [xC,yC] = closest_point(k,Nx.Phi1(k),Ny.Phi1(k),Phi1(k),x,y);
    [xC,yC] = closest_pointInterp(k,Nx.Phi1(k),Ny.Phi1(k),Phi1(k),x,y,Phi1);
    Ghost.Phi1(i).xc=xC;
    Ghost.Phi1(i).yc=yC;
    [stencil,coeffsD,none,nx,ny,coeffs_dx,coeffs_dy] = coeffsLSstencil(x,y,dx,dy,k,Phi1,Mask,xC,yC,1);
    [stencil_c,coeffsD_c,none,none,none,coeffs_dx_c,coeffs_dy_c] = coeffsLSstencil(x,y,dx,dy,k,Phi1,Mask,xC,yC,0);
    Ghost.Phi1(i).G=G(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi1(i).nu=NU(stencil_c(:))'*coeffsD_c(:);
%     Ghost.Phi1(i).Pterm=Pterm(stencil_c(:))'*coeffsD_c(:);
%     Ghost.Phi1(i).Tterm=Tterm(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi1(i).xt1inv=Xt1inv(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi1(i).yt1inv=Yt1inv(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi1(i).stencil=stencil;
    Ghost.Phi1(i).coeffsD=coeffsD;
    Ghost.Phi1(i).coeffs_dx=coeffs_dx;
    Ghost.Phi1(i).coeffs_dy=coeffs_dy;
    Ghost.Phi1(i).stencil_c=stencil_c;
    Ghost.Phi1(i).coeffsD_c=coeffsD_c;
    Ghost.Phi1(i).coeffs_dx_c=coeffs_dx_c;
    Ghost.Phi1(i).coeffs_dy_c=coeffs_dy_c;
    Ghost.Phi1(i).nx=nx;
    Ghost.Phi1(i).ny=ny;
end

for i=1:length(Ghost.Phi2)
    FLAG=2;
    k=Ghost.Phi2(i).index;
%     [xC,yC] = closest_point(k,Nx.Phi2(k),Ny.Phi2(k),Phi2(k),x,y);
    [xC,yC] = closest_pointInterp(k,Nx.Phi2(k),Ny.Phi2(k),Phi2(k),x,y,Phi2);
    Ghost.Phi2(i).xc=xC;
    Ghost.Phi2(i).yc=yC;
    [stencil,coeffsD,none,nx,ny,coeffs_dx,coeffs_dy] = coeffsLSstencil(x,y,dx,dy,k,Phi2,Mask,xC,yC,1);
    [stencil_c,coeffsD_c,none,none,none,coeffs_dx_c,coeffs_dy_c] = coeffsLSstencil(x,y,dx,dy,k,Phi2,Mask,xC,yC,0);
    Ghost.Phi2(i).G=G(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi2(i).nu=NU(stencil_c(:))'*coeffsD_c(:);
%     Ghost.Phi2(i).Pterm=Pterm(stencil_c(:))'*coeffsD_c(:);
%     Ghost.Phi2(i).Tterm=Tterm(stencil_c(:))'*coeffsD_c(:);    
    Ghost.Phi2(i).xt1inv=Xt1inv(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi2(i).yt1inv=Yt1inv(stencil_c(:))'*coeffsD_c(:);    
    Ghost.Phi2(i).stencil=stencil;
    Ghost.Phi2(i).coeffsD=coeffsD;
    Ghost.Phi2(i).coeffs_dx=coeffs_dx;
    Ghost.Phi2(i).coeffs_dy=coeffs_dy;
    Ghost.Phi2(i).stencil_c=stencil_c;
    Ghost.Phi2(i).coeffsD_c=coeffsD_c;
    Ghost.Phi2(i).coeffs_dx_c=coeffs_dx_c;
    Ghost.Phi2(i).coeffs_dy_c=coeffs_dy_c;    
    Ghost.Phi2(i).nx=nx;
    Ghost.Phi2(i).ny=ny;
end

for i=1:length(Ghost.Bdy)
    k=Ghost.Bdy(i).index;
    ig = rem(k-1,nn)+1;
    jg = ceil(k/nn);
    Ghost.Bdy(i).xc=x(jg);
    Ghost.Bdy(i).yc=y(ig);
    type=Ghost.Bdy(i).type;
    ss=((type==0)-(type==1))*nn+((type==2)-(type==3));
    Ghost.Bdy(i).stencil=[k,k+ss,k+2*ss];
    Ghost.Bdy(i).coeffsD=[1,0,0];
    Ghost.Bdy(i).coeffsN=Mask(k+2*ss)*[3 -4 1]/2/dx+(1-Mask(k+2*ss))*[1 -1 0]/dx;
    % STO ASSUMENDO dx==dy!!!
end

