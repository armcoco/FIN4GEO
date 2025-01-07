close all

FontSize=12;
MarkerSize=14;
LineWidth=2.;
MarkerFaceColor=[1 1 1];

subplot(1,1,1,'FontSize',FontSize)
axis square
% axes('Fontsize',Fontsize)
hold on

xx=NN;
yyiS=err;


PiS=polyfit(log10(xx),log10(yyiS),1);

xl=[1.5,4];
yiS=PiS(2)+xl.*PiS(1);

plot(log10(xx),log10(yyiS),'dr','MarkerSize',MarkerSize,'MarkerFaceColor',MarkerFaceColor)
plot(xl,yiS,'-b','LineWidth',LineWidth)
legend('L^\infty error in the solution',['bestfit slope=',num2str(-PiS(1),'%.2f')])
% legend('off')
xlabel('Ln(N)','FontSize',FontSize);
ylabel('Ln(error)','FontSize',FontSize)

% plot2svg('bestfit_flower_1e0_1e6.svg')
