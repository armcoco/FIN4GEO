function vc = coarseFW(vc0,vf,mask)
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

% STENCIL=zeros(nnc*mmc,9); %%aa

% La coarsificazione � del tipo FW (Full-weithing, pag. 43 del Trottemberg)
for jc=2:mmc-1
    for ic=2:nnc-1
        ifine = 2*ic-1;
        jfine = 2*jc-1;
        if mask(ifine,jfine) == 0
            continue;
        end
%         val = 4*vf(ifine,jfine);
        tab = mask(ifine-1:ifine+1,jfine-1:jfine+1);
        
        if tab == ones(3,3);
            val = [1 2 1; 2 4 2; 1 2 1];
        elseif tab(1:2,:) == [1 1 1; 1 1 1];
            val = [ 2 4 2; 2 4 2; 0 0 0 ];
        elseif tab(2:3,:) == [1 1 1 ; 1 1 1];
            val =  [ 0 0 0; 2 4 2; 2 4 2];
        elseif tab(:,1:2) == [1 1 1 ; 1 1 1]';
            val = [ 2 2 0; 4 4 0; 2 2 0 ];
        elseif tab(:,2:3) == [1 1 1 ; 1 1 1]';
            val = [ 0 2 2; 0 4 4; 0 2 2 ];
        elseif tab(1:2,1:2) == [1 1 ; 1 1];
            val = [4 4 0; 4 4 0; 0 0 0];
        elseif tab(1:2,2:3) ==  [1 1 ; 1 1];
            val = [0 4 4 ; 0 4 4; 0 0 0 ];
        elseif tab(2:3,1:2) ==  [1 1 ; 1 1];
            val = [0 0 0; 4 4 0; 4 4 0 ];
        elseif tab(2:3,2:3) ==  [1 1 ; 1 1];
            val = [0 0 0; 0 4 4; 0 4 4 ];
        elseif tab(2,2:3) == [ 1 1];
            val = [0 0 0; 0 8 8; 0 0 0];
        elseif tab(2,1:2) == [ 1 1];
            val = [0 0 0; 8 8 0; 0 0 0];
        elseif tab(1:2,2) == [ 1 1]';
            val = [0 8 0; 0 8 0; 0 0 0];
        elseif tab(2:3,2) == [ 1 1]';
            val = [0 0 0; 0 8 0; 0 8 0];
        else
            val = [0 0 0; 0 16 0; 0 0 0];
        end
        
        coeffs = val.*vf(ifine-1:ifine+1,jfine-1:jfine+1);
        vc(ic,jc) = 1/16*sum(coeffs(:));
        
%         STENCIL(ic+(jc-1)*nnc,:)=val; %%aa
    end
end
        
        
        
        
        
%         if mask(ifine+1,jfine) == 1 && mask(ifine-1,jfine) == 1
%             val = val + vf(ifine+1,jfine) + vf(ifine-1,jfine);
%         else 
%             val = val + 2*vf(ifine,jfine);
%         end
%         if mask(ifine,jfine+1) == 1 && mask(ifine,jfine-1) == 1
%             val = val + vf(ifine,jfine+1) + vf(ifine,jfine-1);
%         else 
%             val = val + 2*vf(ifine,jfine);
%         end
%         vc(ic,jc) = val/8;
%     end
% end

% Jc=2:mmc-1;
% Ic=2:nnc-1;
% Ifine = 2*Ic-1;
% Jfine = 2*Jc-1;
% vc(Ic,Jc) = 1/8*((mask(Ifine,Jfine) == 1).*(...
%     4*vf(Ifine,Jfine) + ...
%     (mask(Ifine+1,Jfine) == 1 & mask(Ifine-1,Jfine) == 1).*(vf(Ifine+1,Jfine) + vf(Ifine-1,Jfine)) + ...
%     (mask(Ifine+1,Jfine) == 0 | mask(Ifine-1,Jfine) == 0).*( 2*vf(Ifine,Jfine) ) + ...
%     (mask(Ifine,Jfine+1) == 1 & mask(Ifine,Jfine-1) == 1).*( vf(Ifine,Jfine+1) + vf(Ifine,Jfine-1) ) + ...
%     (mask(Ifine,Jfine+1) == 0 | mask(Ifine,Jfine-1) == 0).*( 2*vf(Ifine,Jfine) ) ...
%     ));







