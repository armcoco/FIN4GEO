close all

FontSize = 16;
MarkerSize = 8;
LineWidth = 1.;
MarkerFaceColor = [1 1 1];
set(0,'defaultaxesfontsize',FontSize,'defaultlinemarkersize',MarkerSize,...
    'defaultlinelinewidth',LineWidth,'defaultlinemarkerfacecolor',MarkerFaceColor);

figure(1)

subplot(3,1,1)
% % % plot(Xusup,Vtime,'lineWidth',1.1)
% % % xlim([0 10000])
% % % title('Vertical displacement [m]')
% % % xlabel('Radial distance [m]')
% % % ylabel('Vertical displacement [m]')
% % % legend(legend_str)

% perc=Vtime(:,3)./(Vtime(:,2)+Vtime(:,3));
plot(Xusup,perc*100,'lineWidth',1.1)
xlim([0 10000])
ylim([-20 100])
title('Temperature contribution to the vertical displacement [%]')
xlabel('Radial distance [m]')
ylabel('percentage [%]')
legend(legend_str)

% figure(2)
subplot(3,1,2)
solPplot=solP;
% solPplot(Yu>0)=NaN;
% surf(Xu,Yu,solP/1e6)
% shading interp
% view(2)
xu=xt__(x);
yu=yt__(y);
ix = 0-100<=xu & xu<=1e4+100;
iy = -1.5e3-100 <= yu & yu<=0+100;
x_pl=Xu(iy,ix);
y_pl=Yu(iy,ix);
sol_pl=solPplot(iy,ix)/1e6;
contourf(x_pl,y_pl,sol_pl)
% contourf(Xu(iy,ix),Yu(iy,ix),solPplot(iy,ix)/1e6)
axis equal
xlim([0 10000])
ylim([-1500 0])
colorbar
title('pressure change [MPa]')

subplot(3,1,3)
solTplot=solT;
% solTplot(Yu>0)=NaN;
% surf(Xu,Yu,solT)
% shading interp
% view(2)
contourf(Xu(iy,ix),Yu(iy,ix),solTplot(iy,ix))
axis equal
xlim([0 10000])
ylim([-1500 0])
colorbar
title('temperature change [C]')
