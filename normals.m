function [Nx,Ny]=normals(kk,Phi,dx,dy)

[nn,mm]=size(Phi);

Nx = zeros(1,length(kk));
Ny = zeros(1,length(kk));

for pt=1:length(kk)
    xd=2;
    yd=2;
    k=kk(pt);
    if k-nn<1
        Phil = Phi(k);
        xd=1;
    else
        Phil = Phi(k-nn);
    end
    
    if k+nn>nn*mm
        Phir = Phi(k);
        xd=1;
    else
        Phir = Phi(k+nn);
    end
    
    if rem(k,nn)==0
        Phiu = Phi(k);
        yd=1;
    else
        Phiu = Phi(k+1);
    end
    
    if rem(k,nn)==1
        Phid = Phi(k);
        yd=1;
    else
        Phid = Phi(k-1);
    end
    
    Nx(pt) = (Phir-Phil)/(xd*dx);
    Ny(pt) = (Phiu-Phid)/(yd*dy);
end

Nx = reshape(Nx,size(kk));
Ny = reshape(Ny,size(kk));

modulo = sqrt(Nx.^2+Ny.^2);
modulo1 = modulo.*(modulo>0)+(modulo==0);
Nx=Nx./modulo1;
Ny=Ny./modulo1;

Nx=Nx.*(modulo~=0)+(modulo==0);
Ny=Ny.*(modulo~=0);

