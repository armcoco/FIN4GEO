function val=phi2_finite_disl_(X,Y,P1,P2)

[nn,mm]=size(X);
nnmm=nn*mm;
val=zeros(nn,mm);
for k=1:nnmm
    x=X(k);
    y=Y(k);
    P=[x y];
    if OnDislocation(P,P1,P2)
        val(k)=abs(phi2_disl_(x,y,P1,P2));
    else
        d1=sqrt((x-P1(1)).^2+(y-P1(2)).^2);
        d2=sqrt((x-P2(1)).^2+(y-P2(2)).^2);
        val(k)=min(d1,d2);
    end
end
