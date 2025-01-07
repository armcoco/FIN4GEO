function L_SP = L_SP__(x,y)

L_SP = zeros(size(x));

%layer A
k=find(x(:)<=50);
L_SP(k)=-1e-9;

%layer B
k=find(x(:)>50 & x(:)<150);
L_SP(k)=-3.333e-10;

%layer C
k=find(x(:)>150);
L_SP(k)=-1e-10;