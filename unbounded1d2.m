clear all
close all
clc

global c_parameter m_par

iter=1;
NN=[100 200 400 800 1600 3200];
NN=[100 200 400 800];
% NN = ceil([80]*1.2.^[0:20]);
% NN=100;
M=0.1:0.1:2;
M=0.2:0.4:3;
% M=0.2:0.2:1;
% NN=100;
M=1;
ord=[];
for m=M
    err=[];
    for N=NN
        % N=1000;
        
        a=0;
        b=1;
        x=linspace(a,b,N)';
        dx=(b-a)/(N-1);
        
        % au=2*a;
        % bu=2*b;
        % xu=2*x;
        % dxu=(bu-au)/(N-1);
        
        % Au=zeros(N,N);
        % for i=2:N-1
        %     Au(i,[i-1,i,i+1])=[-1 2 -1]/dxu^2;
        % end
        % Au(1,1)=1;
        % Au(N,N)=1;
        % RHSu=-2*ones(N,1); %cos(x);
        % RHSu([1,N])=xu([1,N]).^2;
        %
        % uu=Au\RHSu;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        syms sx;
        sxt=tan(pi/2*sx);
        su_ex=exp(-sx^2);
        sxt=2*sx;
        su_ex=sx^2;
        % sxt=(sx+2)^2;
        % su_ex=sx^2;
        %     sxt=-1e2*log(1-sx);
        %     su_ex=exp(-sx^2);
        %     sxt=sx/(1-sx);
        %     su_ex=exp(-sx^2);
        c=1;
        %     m=1;
        c_parameter=c;
        m_par=m;
        sxt=c*sx/(1-sx^2)^m;
        %     sxt=2*sx;
        %     sxt=sx/(1-sx^2);
        su_ex=exp(-sx^2);
        %     su_ex=1/(sx^2+1);
        su_ex=1/(sx+1); %sx^2+1;
        
%         sxt=sx;
        su_ex=1/(c_parameter*sx/(1-sx^2)^m_par+1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        sxt1=diff(sxt,1);
        sxt2=diff(sxt,2);
        sf=-diff(su_ex,2);
        sf=-1/(sxt1)*diff(1/sxt1*diff(su_ex,1),1);
        sxt=sx;
        matlabFunction(sxt,'file','xt_.m','vars',[sx]);
        %     matlabFunction(sxt1,'file','xt1_.m','vars',[sx]);
        %     matlabFunction(sxt2,'file','xt2_.m','vars',[sx]);
        matlabFunction(su_ex,'file','u_ex_.m','vars',[sx]);
        matlabFunction(sf,'file','f_.m','vars',[sx]);
        sF=f_(sxt);
        matlabFunction(sF,'file','F_.m','vars',[sx]);
        
        xu=xt_(x);
        u_exu=u_ex_(xu);
        
        A=zeros(N,N);
        for i=2:N-1
            %         A(i,[i-1,i,i+1])=[-1 2 -1]/dx^2+xt2_(x(i))/xt1_(x(i))*[-1 0 1]/2/dx;
            gammaL=1/xt1_((x(i-1)+x(i))/2);
            gammaR=1/xt1_((x(i+1)+x(i))/2);
            A(i,[i-1,i,i+1])=1/xt1_(x(i))*[-gammaL gammaL+gammaR -gammaR]/dx^2;
            %         gammaL=dx/(xt_(x(i))-xt_(x(i-1)));
            %         gammaR=dx/(xt_(x(i+1))-xt_(x(i)));
            %         A(i,[i-1,i,i+1])=2*dx/(xt_(x(i+1))-xt_(x(i-1)))*[-gammaL gammaL+gammaR -gammaR]/dx^2;
        end
        A(1,1)=1;
        A(N,N)=1;
        % RHS=-2*4*ones(N,1);
        RHS=zeros(N,1);
        %     RHS(:)=F_(x); %.*xt1_(x).^2;
        RHS(:)=f_(xu); %.*xt1_(x).^2;
        RHS([1,N])=u_exu([1,N]);
        
        u=A\RHS;
        
        df=u(2:end-1)-u_exu(2:end-1);
        err(end+1)=max(abs(df(:)))
            p_norm=2;
%             err(end+1)=norm(df,p_norm)/norm(u,p_norm)
        %     err(end+1)=(sum(abs(df).^p_norm.*xt1_(x(2:end-1)))./sum(abs(u(2:end-1)).^p_norm.*xt1_(x(2:end-1)))).^(1/p_norm)
%         iter=iter+1;
    end
    bestfit1d;
    ord(end+1)=PiS(1)
end





