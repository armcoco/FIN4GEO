function [Ghost]=setGrid_P_and_T(Ghost,Pterm,Tterm)

for i=1:length(Ghost.Phi1)
    stencil_c=Ghost.Phi1(i).stencil_c;
    coeffsD_c=Ghost.Phi1(i).coeffsD_c;
    Ghost.Phi1(i).Pterm=Pterm(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi1(i).Tterm=Tterm(stencil_c(:))'*coeffsD_c(:);
end

for i=1:length(Ghost.Phi2)
    stencil_c=Ghost.Phi2(i).stencil_c;
    coeffsD_c=Ghost.Phi2(i).coeffsD_c;
    Ghost.Phi2(i).Pterm=Pterm(stencil_c(:))'*coeffsD_c(:);
    Ghost.Phi2(i).Tterm=Tterm(stencil_c(:))'*coeffsD_c(:);
end

