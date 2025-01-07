function Sig = sig__(x,y)

Sig = zeros(size(x));

%layer A
k=find(x(:)<=50);
Sig(k)=1;

%layer B
k=find(x(:)>50 & x(:)<150);
Sig(k)=0.3333;

%layer C
k=find(x(:)>150);
Sig(k)=0.1;


