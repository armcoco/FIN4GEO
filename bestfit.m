close all

xx=N;
yy=LL1;

P1=polyfit(log10(xx),log10(yy),1);

loglog(xx,yy,'*r')
hold on
xl=[10,1000];
yl=xl.^P1(1)*10.^P1(2);
plot(xl,yl,'b')
legend('off')

% plot2svg('C:\Users\Armando\Documents\Dottorato\articoli\articolo\articolo\Figures\bestfitFFPM1.svg')


figure

xx=N;
yy=LLi;

Pi=polyfit(log10(xx),log10(yy),1);

loglog(xx,yy,'*r')
hold on
xl=[10,1000];
yl=xl.^Pi(1)*10.^Pi(2);
plot(xl,yl,'b')
legend('off')

% plot2svg('C:\Users\Armando\Documents\Dottorato\articoli\articolo\articolo
% \Figures\bestfitFFPMi.svg')