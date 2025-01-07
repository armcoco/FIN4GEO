function [stencil,coeffsD,coeffsN,nx,ny,coeffs_ux,coeffs_uy] = coeffsLSstencil(x,y,dx,dy,k,Phi,Mask,xC,yC,upwind)
%
% Ricostruisce la soluzione numerica nel punto (xC,yC) utilizzando uno
% stencil contenente l'indice k e formato da s*s punti in direzione Upwind
%

global FLAG;

s=3;
nn=length(y);
mm=length(x);

ii=0; 
jj=0;
ig = rem(k-1,nn)+1;
jg = ceil(k/nn);
xG=x(jg); yG=y(ig);
sx=sign(xC-xG);  sy=sign(yC-yG);
if sx==0
    sx=sgn(Phi(k-nn)-Phi(k+nn));
    sx=-1;
    sx=1*(jg<=mm/2)-1*(jg>mm/2);
end
if sy==0
    sy=sgn(Phi(k-1)-Phi(k+1));
    sy=-1;
    sy=1*(ig<=nn/2)-1*(ig>nn/2);
end

stencil=zeros(3,3);
% epsl=1e-1;
% modl=sqrt((xC-xG)^2+(yC-yG)^2);
if (jg==2 && sx<0) || (jg==mm-1 && sx>0) || upwind==0 %abs(xC-xG)/modl<epsl
    jj=1;
end
if (ig==2 && sy<0) || (ig==nn-1 && sy>0) || upwind==0 %abs(yC-yG)/modl<epsl
    ii=1;
end

for j=1:3
    for i=1:3
        stencil(i+3*(j-1))=k+sx*(j-1-jj)*nn+sy*(i-1-ii);
    end
end
    
% stencil=k+sx*nn.*repmat(0:s-1,s,1)+sy.*repmat(0:s-1,s,1)';

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

% TOLL=0.1;
% if abs(thty-1/2)<TOLL && thtx<TOLL
%     weights_y=[1-thty 0 thty];
% end
% if abs(thtx-1/2)<TOLL && thty<TOLL
%     weights_x=[1-thtx 0 thtx];
% end

coeffsD=weights_y'*weights_x;
coeffs_dx=sx*weights_y'*weights_dx;
coeffs_dy=sy*weights_dy'*weights_x;
nx=sum(sum(Phi(stencil).*coeffs_dx));
ny=sum(sum(Phi(stencil).*coeffs_dy));
module=sqrt(nx^2+ny^2);
nx=nx/module;
ny=ny/module;
% if FLAG==1
%     [nx,ny] = normal_phi1_(xt_(xC),yt_(yC));
% else
%     [nx,ny] = normal_phi2_(xt_(xC),yt_(yC));
% end

%%%%%%%%%% check whether or not I have to shift the stencil
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
    %     disp(['riduco: ',num2str(k)])
    %                     pause
    stencil=k+sx*nn.*repmat(0:s-1,s,1)+sy.*repmat(0:s-1,s,1)';
    coeffsD=[1 0 0;0 0 0; 0 0 0];
    coeffs_ux=sx*([-1 1 0;0 0 0; 0 0 0]/dx*(Mask(stencil(4))==1)+[0 0 0; -1 1 0; 0 0 0]/dx*(Mask(stencil(4))==0));
    coeffs_uy=sy*([-1 0 0; 1 0 0; 0 0 0]/dy*(Mask(stencil(2))==1)+[0 -1 0; 0 1 0; 0 0 0]/dy*(Mask(stencil(2))==0));
    coeffsN=coeffs_ux*nx+coeffs_uy*ny;
    return;
end

%%%%%%%%%%%%%%%%%%%%%

coeffs_ux=coeffs_dx;
coeffs_uy=coeffs_dy;
coeffsN=coeffs_ux*nx+coeffs_uy*ny;

