close all
figure
hold on
axis equal
a=Xu(nn+1); b=Xu(nn*mm-nn); c=Yu(2); d=Yu(nn-1);
xlim([a b])
ylim([c d])
% plot(X,Y,'.r')

plot([Xu(1),Xu(1)],[Yu(1),Yu(end)],'r','lineWidth',2)
plot([Xu(1),Xu(end)],[Yu(end),Yu(end)],'r','lineWidth',2)
plot([Xu(end),Xu(end)],[Yu(end),Yu(1)],'r','lineWidth',2)
plot([Xu(end),Xu(1)],[Yu(1),Yu(1)],'r','lineWidth',2)

contour(Xu,Yu,Phi2,[0 0],'r','lineWidth',2)

plot(Xu(Inside.all),Yu(Inside.all),'g.')
% % % plot(Xu([Ghost.Phi2.index]),Yu([Ghost.Phi2.index]),'ok')
plot(Xu([Ghost.Bdy.index]),Yu([Ghost.Bdy.index]),'or')

if ~DomInf
    contour(Xu,Yu,Phi1,[0 0],'r','lineWidth',2)
    for i=1:length(Ghost.Phi1)
        k=Ghost.Phi1(i).index;
        plot([Xu(k),xt__(Ghost.Phi1(i).xc)],[Yu(k),yt__(Ghost.Phi1(i).yc)],'b')
        plot(xt__(Ghost.Phi1(i).xc),yt__(Ghost.Phi1(i).yc),'.r')
    end
    plot(Xu([Ghost.Phi1.index]),Yu([Ghost.Phi1.index]),'ob')
    quiver(xt__([Ghost.Phi1.xc]),yt__([Ghost.Phi1.yc]),[Ghost.Phi1.nx]*dx,[Ghost.Phi1.ny]*dx,0)
end

% % % for i=1:length(Ghost.Phi2)
% % %     k=Ghost.Phi2(i).index;
% % %     plot([Xu(k),xt__(Ghost.Phi2(i).xc)],[Yu(k),yt__(Ghost.Phi2(i).yc)],'b')
% % %     plot(xt__(Ghost.Phi2(i).xc),yt__(Ghost.Phi2(i).yc),'.r')
% % % end


% sol_i = sol_i_(X,Y); % exact solution in matrix fashion
% sol_e = sol_e_(X,Y); % exact solution in matrix fashion
% [sol_i,sol_e]=controlDim(nn,mm,sol_i,sol_e);
% ErrorArray=abs([u_i(Inside.i)-sol_i(Inside.i) ; u_e(Inside.e)-sol_e(Inside.e)]);
% ErrorMatrix=zeros(nn,mm);
% ErrorMatrix([Inside.i;Inside.e])=ErrorArray;
% contour(X,Y,ErrorMatrix);
