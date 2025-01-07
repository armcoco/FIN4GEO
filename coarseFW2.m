function vc = coarseFW2(vc0,vf,mask,Nx,Ny)
%
% E' l'operatore di transfer da una griglia fine ad una rada, settato
% soprattutto per il residuo nei punti interni.
% vc0 � il vettore nella griglia rada iniziale (tale input serve quando ad
% esempio si deve coarsificare il residuo sia nei punti interni, che in
% quelli esterni per le condizioni di Dirichlet: in questo modo non si
% perde la precedente coarsificazione).
% vf � il vettore nella griglia fine
% mask � una matrice i cui elementi sono 1 (se il punto � interno alla zona
% che si vuole coarsificare) e 0 (se � esterno).
%

vc=vc0;
[nn,mm] = size(vf);
n=nn-1;
m=mm-1;
nc=n/2;
mc=m/2;
nnc=nc+1;
mmc=mc+1;

% La coarsificazione � del tipo FW (Full-weithing, pag. 43 del Trottemberg)
for jc=2:mmc-1
    for ic=2:nnc-1
        ifine = 2*ic-1;
        jfine = 2*jc-1;
        k=ifine+(jfine-1)*nn;
        if mask(ifine,jfine) == 0
            continue;
        end
        nx=Nx(k);
        ny=Ny(k);
        sx=sgn(nx); sy=sgn(ny);
        si1=sy*(abs(nx)>abs(ny))+sx*nn*(abs(nx)<=abs(ny));
        si2=sx*nn*(abs(nx)>abs(ny))+sy*(abs(nx)<=abs(ny));
        stencil=zeros(3,3);
        for i2=-1:1 
            for i1=-1:1
                stencil((i1+1)+3*(i2+1)+1)=k+si1*i1+si2*i2;
            end
        end
        val=vf(stencil);
        tab = mask(stencil);
%         if ~all(tab(:)) | 1==1
%                 val=vf(k)*ones(3,3);
%         end 
        if ~all(tab([7,8,9]))
            val([7,8,9])=val([1,2,3]);
        end
        if ~all(tab([3,6]))
            val([3,6,9])=val([1,4,7]);
        end
        if ~all(tab([1,4]))
            val([1,4,7])=val([2,5,8]);
            val([3,6,9])=val([2,5,8]);
        end
        if ~all(tab([2]))
            val([1,2,3])=val([4,5,6]);
            val([7,8,9])=val([4,5,6]);
        end
        coeffs = [1 2 1; 2 4 2; 1 2 1]/16;
        vc(ic,jc) = sum(coeffs(:).*val(:));
%         vc(ic,jc) = vf(k);
    end
end
