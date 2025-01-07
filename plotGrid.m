close all
figure
hold on
axis equal
a=X(1); b=X(nn*mm); c=Y(1); d=Y(nn);
xlim([a b])
ylim([c d])
% plot(X,Y,'.r')

plot([X(1),X(1)],[Y(1),Y(end)],'r','lineWidth',2)
plot([X(1),X(end)],[Y(end),Y(end)],'r','lineWidth',2)
plot([X(end),X(end)],[Y(end),Y(1)],'r','lineWidth',2)
plot([X(end),X(1)],[Y(1),Y(1)],'r','lineWidth',2)

contour(X,Y,Phi1,[0 0],'r','lineWidth',2)
contour(X,Y,Phi2,[0 0],'r','lineWidth',2)

plot(X(Inside.all),Y(Inside.all),'g.')
plot(X([Ghost.Phi1.index]),Y([Ghost.Phi1.index]),'ob')
plot(X([Ghost.Phi2.index]),Y([Ghost.Phi2.index]),'ok')
plot(X([Ghost.Bdy.index]),Y([Ghost.Bdy.index]),'or')

for i=1:length(Ghost.Phi1)
    k=Ghost.Phi1(i).index;
    plot([X(k),Ghost.Phi1(i).xc],[Y(k),Ghost.Phi1(i).yc],'b')
    plot(Ghost.Phi1(i).xc,Ghost.Phi1(i).yc,'.r')
end

for i=1:length(Ghost.Phi2)
    k=Ghost.Phi2(i).index;
    plot([X(k),Ghost.Phi2(i).xc],[Y(k),Ghost.Phi2(i).yc],'b')
    plot(Ghost.Phi2(i).xc,Ghost.Phi2(i).yc,'.r')
end

quiver([Ghost.Phi1.xc],[Ghost.Phi1.yc],[Ghost.Phi1.nx]*dx,[Ghost.Phi1.ny]*dx,0)

% sol_i = sol_i_(X,Y); % exact solution in matrix fashion
% sol_e = sol_e_(X,Y); % exact solution in matrix fashion
% [sol_i,sol_e]=controlDim(nn,mm,sol_i,sol_e);
% ErrorArray=abs([u_i(Inside.i)-sol_i(Inside.i) ; u_e(Inside.e)-sol_e(Inside.e)]);
% ErrorMatrix=zeros(nn,mm);
% ErrorMatrix([Inside.i;Inside.e])=ErrorArray;
% contour(X,Y,ErrorMatrix);
