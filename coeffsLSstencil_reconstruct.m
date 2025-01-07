function [stencil,coeffsD,coeffsN,coeffs_ux,coeffs_uy,nx,ny] = coeffsLSstencil_reconstruct(x,y,dx,dy,k,Phi,Mask,xC,yC)
%
% Ricostruisce la soluzione numerica nel punto (xB,yB) utilizzando uno
% stencil che racchiude (xB,yB) e formato da s*s punti in direzione Upwind
%

nn=length(y);
mm=length(x);

ig = rem(k-1,nn)+1;
jg = ceil(k/nn);
xG=x(jg); yG=y(ig);
sx=sign(xG-xC);  sy=sign(yG-yC); %%% attenzione! Il segno è cambiato rispetto a coeffsLSstencil
if sx==0
    sx=1; %sgn(Phi(k-nn)-Phi(k+nn));
end
if sy==0
    sy=1; %sgn(Phi(k-1)-Phi(k+1));
end

stencil=k+sx*nn.*repmat(-1:1,3,1)+sy*repmat(-1:1,3,1)';

k1=stencil(1,1);
ig = rem(k1-1,nn)+1;
jg = ceil(k1/nn);
xG1=x(jg); yG1=y(ig);
thtx=abs(xC-xG1)/2/dx;
thty=abs(yC-yG1)/2/dy;

weights_x=[(1-2*thtx)*(1-thtx) 4*thtx*(1-thtx) thtx*(2*thtx-1)];
weights_y=[(1-2*thty)*(1-thty) 4*thty*(1-thty) thty*(2*thty-1)];
weights_dx=( [-1 1 0] + [1 -2 1]*(2*thtx-0.5) )/dx;
weights_dy=( [-1 1 0] + [1 -2 1]*(2*thty-0.5) )/dy;

coeffsD=weights_y'*weights_x;
coeffs_dx=sx*weights_y'*weights_dx;
coeffs_dy=sy*weights_dy'*weights_x;

nx=sum(sum(Phi(stencil).*coeffs_dx));
ny=sum(sum(Phi(stencil).*coeffs_dy));
module=sqrt(nx^2+ny^2);
nx=nx/module;
ny=ny/module;
[nx,ny] = normal_phi1_(xC,yC);

%%%%%%%%%% check whether or not I have to slide the stencil
if abs(xC-xG)>abs(yC-yG)
    direction=1;
    s1=1;
else
    direction=0;
    s1=3;
end
    
sxy=sx*(direction==1)+sy*(direction==0);

for i=0:2
    Face(i+1)=1+i*s1;
end

for ind_face=Face
    for pt=0:1
        ind=ind_face+pt*3^direction;
        if Mask(stencil(ind))
            break;
        else
            stencil(ind)=stencil(ind)+sxy*3*nn^direction;
            cD=coeffsD(ind);    cdx=coeffs_dx(ind);    cdy=coeffs_dy(ind);
            coeffsD(ind)=0;     coeffs_dx(ind)=0;      coeffs_dy(ind)=0;
            c_extrap=[1 3 -3]*(pt==0)+[-3 1 3]*(pt==1);
            coeffsD(ind_face+(0:2)*3^direction)=coeffsD(ind_face+(0:2)*3^direction)+cD*c_extrap;
            coeffs_dx(ind_face+(0:2)*3^direction)=coeffs_dx(ind_face+(0:2)*3^direction)+cdx*c_extrap;
            coeffs_dy(ind_face+(0:2)*3^direction)=coeffs_dy(ind_face+(0:2)*3^direction)+cdy*c_extrap;
        end
    end
end

if ~all(Mask(stencil(:)))
%     disp(['riduco_r: ',num2str(k)])
    %                     pause
    stencil=k+sx*nn.*repmat(-1:1,3,1)+sy*repmat(-1:1,3,1)';
    coeffsD=[0 0 0;0 1 0; 0 0 0];
    coeffs_ux=sx*[0 0 0;-1 1 0; 0 0 0]/dx;
    coeffs_uy=sy*[0 -1 0;0 1 0; 0 0 0]/dy; 
    coeffsN=coeffs_ux*nx+coeffs_uy*ny;
    return;
end

%%%%%%%%%%%%%%%%%%%%%

coeffs_ux=coeffs_dx;
coeffs_uy=coeffs_dy;
coeffsN=coeffs_ux*nx+coeffs_uy*ny;
