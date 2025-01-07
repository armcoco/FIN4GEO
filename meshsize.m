xL=0;
xR=1e4;
zL=-1500;
zR=0;

Xcenter=X(1:l);
Zcenter=Z(1:l:end);

xsize=zeros(l,1);
xsize(1)=(Xcenter(1)-xL)*2;
for i=2:l
    xsize(i)=(Xcenter(i)-(Xcenter(i-1)+xsize(i-1)/2))*2;
end

zsize=zeros(m,1);
zsize(1)=(Zcenter(1)-zL)*2;
for i=2:m
    zsize(i)=(Zcenter(i)-(Zcenter(i-1)+zsize(i-1)/2))*2;
end

[Xsize,Zsize]=meshgrid(xsize,zsize);
Xsize(1,:)=[];
Zsize(:,1)=[];


%%%%%%

xu=zeros(1,l-1);
xu(1:l-1)=Xcenter(1:l-1)+xsize(1:l-1)/2;
zu=Zcenter;
[X_u,Z_u]=meshgrid(xu,zu);
xv=Xcenter;
zv=zeros(1,m-1);
zv(1:m-1)=Zcenter(1:m-1)+zsize(1:m-1)/2;
[X_v,Z_v]=meshgrid(xv,zv);
