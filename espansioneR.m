function res = espansioneR(res,expres,nx,ny,direction)

[nn,mm] = size(res);
% n=nn-1;
% m=mm-1;
% inside=find(mask);

k=expres;
% [nx,ny]=gradPhi1_(X(expres),Y(expres));
modulo=1; %nx.^2+ny.^2;
nx=(direction*nx./modulo);
ny=(direction*ny./modulo);
% quadrante = 1*(nx<=0 & ny<0) + 2*(nx>0 & ny<0) + 3*(nx>0 & ny>=0) + 4*(nx<=0 & ny>=0);
% rx = (res(k)-res(k-nn)).*(quadrante==2 | quadrante==3) + (res(k+nn)-res(k)).*(quadrante==1 | quadrante==4);
% ry = (res(k)-res(k-1)).*(quadrante==3 | quadrante==4) + (res(k+1)-res(k)).*(quadrante==1 | quadrante==2);
%     rx = 0.5*(3*res(k)-4*res(k-nn)+res(k-2*nn)).*(quadrante==2 | quadrante==3) - 0.5*(3*res(k)-4*res(k+nn)+res(k+2*nn)).*(quadrante==1 | quadrante==4);
%     ry = 0.5*(3*res(k)-4*res(k-1)+res(k-2)).*(quadrante==3 |
%     quadrante==4) - 0.5*(3*res(k)-4*res(k+1)+res(k+2)).*(quadrante==1 |
%     quadrante==2);

for i=1:length(expres)
    k=expres(i);
    rx = (res(k)-res(k-nn)).*(nx(k)>0) + (res(k+nn)-res(k)).*(nx(k)<=0);
    ry = (res(k)-res(k-1)).*(ny(k)>0) + (res(k+1)-res(k)).*(ny(k)<=0);
    res(k) = res(k) - 0.9/sqrt(2)*(rx.*nx(k)+ry.*ny(k));
end
