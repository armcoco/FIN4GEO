close all

FontSize=12;
MarkerSize=14;
LineWidth=2.;
MarkerFaceColor=[1 1 1];

subplot(1,1,1,'FontSize',FontSize)
axis square
% axes('Fontsize',Fontsize)
hold on

xx=N;
yyiG=LLiG;
yy1G=LL1G;
yyiS=LLiS;
yy1S=LL1S;

yyiG=LL1Su;
yy1G=LL1Sv;
yyiS=LLiSu;
yy1S=LLiSv;

PiG=polyfit(log10(xx),log10(yyiG),1);
PiS=polyfit(log10(xx),log10(yyiS),1);
P2G=polyfit(log10(xx),log10(yy1G),1);
P2S=polyfit(log10(xx),log10(yy1S),1);

xl=[1,3];
% xlim(xl);
% ylim([-15 -1]);
yiG=PiG(2)+xl.*PiG(1);
yiS=PiS(2)+xl.*PiS(1);
y2G=P2G(2)+xl.*P2G(1);
y2S=P2S(2)+xl.*P2S(1);

plot(log10(xx),log10(yyiG),'*r','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor)
plot(xl,yiG,'-.b','LineWidth',LineWidth)
plot(log10(xx),log10(yy1G),'dr','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor)
plot(xl,y2G,'-b','LineWidth',LineWidth)
plot(log10(xx),log10(yyiS),'or','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor)
plot(xl,yiS,'--b','LineWidth',LineWidth)
plot(log10(xx),log10(yy1S),'+r','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor)
plot(xl,y2S,':b','LineWidth',LineWidth)
legend('L^\infty error in the gradient',['bestfit slope=',num2str(-PiG(1),'%.2f')],...
    'L^1 error in the gradient',['bestfit slope=',num2str(-P2G(1),'%.2f')],...
    'L^\infty error in the solution',['bestfit slope=',num2str(-PiS(1),'%.2f')],...
    'L^1 error in the solution',['bestfit slope=',num2str(-P2S(1),'%.2f')])
% legend('off')
xlabel('Ln(N)','FontSize',FontSize);
ylabel('Ln(error)','FontSize',FontSize)

% plot2svg('bestfit_flower_1e0_1e6.svg')
