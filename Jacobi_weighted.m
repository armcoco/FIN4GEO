clear all
close all

N=100;
omega=0.5;
x=linspace(-1,1,N);
u=2*rand(1,N)-1;
u(1)=0;
u(end)=0;


for iter=1:1000
    uold=u;
    plot(x,u,'-*b');
    for i=2:N-1
        u(i)=(uold(i-1)+uold(i+1))/2;
    end
    u=omega*u+(1-omega)*uold;
    
    % plot
    ylim([-1 1])
    text(0,0.7,['\color{black} \fontsize{36} iteration #',num2str(iter)]);
    hold on
    plot([-1 1],[0 0],'r')
    drawnow
    cla
end

