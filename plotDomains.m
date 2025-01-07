close all

figure
hold on
axis square
box on
contourf(X,Y,-Phi2,[0 0],'b','lineWidth',3)
contourf(X,Y,-Phi1,[0 0],'r','lineWidth',3)
colormap([0 1 1])