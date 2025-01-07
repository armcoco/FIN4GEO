clear all
close all
x=linspace(-1,1,1000);
y=x;
[X,Y]=meshgrid(x,y);
phi=X.^2+Y.^2-0.5^2;
plot3(X,Y,phi)

N=100;
x=linspace(-1,1,N);
u0=1*(2*rand(1,N)-1);
load u0.txt;
u=u0;
% u=rand(1,N);
u(1)=0;
u(end)=0;
mm=max(abs(u));

writerObj = VideoWriter('GaussSeidel.avi');
writerObj.FrameRate = 2;
open(writerObj);

% f=5*exp(-x.^2/1e-2);
f=0*(5*exp(-(x-0.5).^2/1e-2)+5*exp(-(x+0.5).^2/1e-2));
% f(-0.5<x & x<0.5)=0.5;
dx=2/1000;

for iter=1:50
    uold=u;
    plot(x,u,'-*b');
    for i=2:N-1
        u(i)=(u(i-1)+u(i+1)+f(i)*dx^2)/2;
    end
%     axis square
%     axis equal
%     xlim([-1 1])
    mm=max(abs(u));
    ylim([-1 1])
    text(0,0.7,['\color{black} \fontsize{36} iteration #',num2str(iter)]);
    hold on
    plot([-1 1],[0 0],'r')
    drawnow
%     pause
    frame = getframe;
   writeVideo(writerObj,frame);
   cla
%     M(iter) = getframe(gca);
end

close(writerObj);
