function [gp,M,RHS,GP]=LinearSystemSparseEllipticAS(dx,dy,RR,Inside,Ghost,F,Xt1inv,Yt1inv)
global DomInf
% w=1; %4/5; %peso dell'iterazione di Jacobi
[nn,mm]=size(F.gp);

% M=speye(2*nn*mm);
RHS=nan*zeros(nn*mm,1);
% Pinv_diag=zeros(2*nn*mm,1);
nnmm=nn*mm;

auLR=Xt1inv;
auBT=Yt1inv;

Outside=[1 nn nn*mm-nn+1 nn*mm];
NNZ=numel(Outside)+5*(nn-2)*(mm-2)+3*(2*(nn-2)+2*(mm-2));
ROWS=zeros(1,NNZ);
COLS=zeros(1,NNZ);
VALS=zeros(1,NNZ);
ROWS(1:numel(Outside))=[Outside];
COLS(1:numel(Outside))=[Outside];
VALS(1:numel(Outside))=1;
cont=numel(Outside);
for i=2:nn-1
    for j=2:mm-1
        k=i+(j-1)*nn;
        gamma1=Xt1inv(k)/RR.centre(k);
        gamma2=Yt1inv(k);
        auR=0.5*gamma1*(auLR(k+nn)+auLR(k))*RR.right(k);
        auL=0.5*gamma1*(auLR(k-nn)+auLR(k))*RR.left(k);
        auT=0.5*gamma2*(auBT(k+1)+auBT(k));
        auB=0.5*gamma2*(auBT(k-1)+auBT(k));
        ROWS(cont+(1:5))=k;
        COLS(cont+(1:5))=[k k-nn k+nn k-1 k+1];
        VALS(cont+(1:5))=[auR+auL+auT+auB -auL -auR -auB -auT]/dx^2;
        cont=cont+5;
        RHS(k)=F.gp(k);
        %     Pinv_diag(k)=dx^2/(auR+auL+auT+auB+a0);
    end
end

for j=2:mm-1
    i=1;
    k = i+(j-1)*nn;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[k k+1 k+2];
    VALS(cont+(1:3))=[1 0 0];
    cont=cont+3;
    RHS(k)=F.gp(k);
    %     Pinv_diag(k)=0.8;
    i=nn;
    k = i+(j-1)*nn;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[k k-1 k-2];
    VALS(cont+(1:3))=[1 0 0];
    cont=cont+3;
    RHS(k)=F.gp(k);
    %     Pinv_diag(k)=0.8;    
end

for i=2:nn-1
    j=1;
    k = i+(j-1)*nn;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[k k+nn k+2*nn];
    VALS(cont+(1:3))=Xt1inv(k)*[3 -4 1]/2/dx;
    cont=cont+3;
    RHS(k)=F.gp(k);
    %     Pinv_diag(k)=0.8;
    j=mm;
    k = i+(j-1)*nn;
    ROWS(cont+(1:3))=k;
    COLS(cont+(1:3))=[k k-nn k-2*nn];
    VALS(cont+(1:3))=[1 0 0];
    cont=cont+3;
    RHS(k)=F.gp(k);
    %     Pinv_diag(k)=0.8;
end

M = sparse(ROWS,COLS,VALS);

GP=M\RHS;
gp=reshape(GP,nn,mm);

% Pinv=diag(Pinv_diag);
% B=speye(2*nnmm)-Pinv*M;
