function vc = coarseE(vf)

[nn,mm] = size(vf);
n=nn-1;
m=mm-1;
nc=n/2;
mc=m/2;
nnc=nc+1;
mmc=mc+1;
vc=zeros(nnc,mmc);

for jc=1:mmc
    for ic=1:nnc
        ifine = 2*ic-1;
        jfine = 2*jc-1;
        vc(ic,jc) = vf(ifine,jfine);
    end
end

% Inc = 2:mmc-1;
% cc = fix(sum(mask(:,2))/2)+2;
% ff = sum(mask(:,2))+2;
% vc(cc,Inc) = (vf(ff,2*Inc-2) + 2*vf(ff,2*Inc-1) + vf(ff,2*Inc) )/4;
% 
% % vc(cc,[2,mmc-1]) = vf(ff,[2,mm-1]); 
% 
% % Inc = 2:mmc-1;
% % vc(6,Inc) = vf(11,2*Inc-1); 