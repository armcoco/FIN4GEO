function [xC,yC,phi_x,phi_y,modulo] = closest_point(kk,nx,ny,R,x,y)
%
% Trova il punto pi� vicino dell'interfaccia del punto k-esimo.
% In ingresso, n � il numero di intervalli in cui divido l'intervallo
% sull'asse y [c,d].

% i = ceil((k-1)/n);
% j = k-n*(i-1);
% x0 = x(i);
% y0 = y(j);

nn=length(y);
mm=length(x);

dx = x(2) - x(1);
dy = y(2) - y(1);

ig = rem(kk-1,nn)+1;
jg = ceil(kk/nn);
xP = x(jg);
yP = y(ig);
xC = xP' - R.*nx;
yC = yP' - R.*ny;


% % % function [xC,yC,quadrante, st] = closest_point(X,Y,Phi,k,n,delta_x,delta_y)
% % % %
% % % % Trova il punto pi� vicino dell'interfaccia del punto k-esimo.
% % % % In ingresso, n � il numero di intervalli in cui divido l'intervallo
% % % % sull'asse y [c,d].
% % % 
% % % ENO = 0;
% % % INTERP = 0;
% % % 
% % % phi_x = zeros(size(X));
% % % phi_y=phi_x;
% % % st = zeros(length(k),6);
% % % 
% % % % i = ceil((k-1)/n);
% % % % j = k-n*(i-1);
% % % % x0 = x(i);
% % % % y0 = y(j);
% % % 
% % % % % Commentare e scommentare se si vuole una ricostruzione ENO per le derivate phi_x e
% % % % phi_y e commentare le righe di codice subito dopo;
% % % if ENO==1 
% % % %     st = zeros(length(k),6);
% % %     kk=k.*(abs(Phi(k)-Phi(k-n))<abs(Phi(k)-Phi(k+n)))+(k-n).*(abs(Phi(k)-Phi(k-n))>=abs(Phi(k)-Phi(k+n)));
% % %     segno=1.*(abs(Phi(k)-Phi(k-n))<abs(Phi(k)-Phi(k+n)))-1.*(abs(Phi(k)-Phi(k-n))>=abs(Phi(k)-Phi(k+n)));
% % %     Q1=(Phi(kk+n)-Phi(kk))/delta_x;
% % %     
% % %     D2kk=0.5/delta_x^2*(Phi(kk+n)-2*Phi(kk)+Phi(kk-n));
% % %     D2kk1=0.5/delta_x^2*(Phi(kk+2*n)-2*Phi(kk+n)+Phi(kk));
% % %     
% % %     cc=D2kk.*(abs(D2kk)<=abs(D2kk1))+D2kk1.*(abs(D2kk)>abs(D2kk1));
% % %     for i=1:length(D2kk)
% % %         if abs(D2kk(i))<=abs(D2kk1(i))
% % %             st(i,1:3)=[kk(i)+n, kk(i), kk(i)-n];
% % %         else
% % %             st(i,1:3)=[ kk(i)+2*n, kk(i)+n, kk(i)];
% % %         end
% % %     end
% % %     Q2=segno.*cc.*delta_x;
% % %     
% % %     phi_x = Q1+Q2;
% % %     
% % %     kk=k.*(abs(Phi(k)-Phi(k-1))<abs(Phi(k)-Phi(k+1)))+(k-1).*(abs(Phi(k)-Phi(k-1))>=abs(Phi(k)-Phi(k+1)));
% % %     segno=1.*(abs(Phi(k)-Phi(k-1))<abs(Phi(k)-Phi(k+1)))-1.*(abs(Phi(k)-Phi(k-1))>=abs(Phi(k)-Phi(k+1)));
% % %     Q1=(Phi(kk+1)-Phi(kk))/delta_y;
% % %     
% % %     D2kk=0.5/delta_y^2*(Phi(kk+1)-2*Phi(kk)+Phi(kk-1));
% % %     D2kk1=0.5/delta_y^2*(Phi(kk+2)-2*Phi(kk+1)+Phi(kk));
% % %     
% % %     cc=D2kk.*(abs(D2kk)<=abs(D2kk1))+D2kk1.*(abs(D2kk)>abs(D2kk1));
% % %     for i=1:length(D2kk)
% % %         if abs(D2kk(i))<=abs(D2kk1(i))
% % %             st(i,4:6)=[kk(i)+1, kk(i), kk(i)-1];
% % %         else
% % %             st(i,4:6)=[ kk(i)+2, kk(i)+1, kk(i)];
% % %         end
% % %     end
% % %     Q2=segno.*cc.*delta_y;
% % %     
% % %     phi_y = Q1+Q2;
% % % else
% % %     phi_x = (Phi(k+n)-Phi(k-n))/(2*delta_x);
% % %     phi_y = (Phi(k+1)-Phi(k-1))/(2*delta_y);
% % % end
% % % 
% % % if INTERP==1
% % %     xC=zeros(size(k));
% % %     yC=zeros(size(k));
% % %     for ii=1:length(k)
% % %         ki=k(ii);
% % %         grigliax = (ki-2*n:n:ki).*(phi_x(ii)>0) + (ki:n:ki+2*n).*(phi_x(ii)<=0);
% % %         griglia = [grigliax-2; grigliax-1; grigliax].*(phi_y(ii)>0) + [grigliax; grigliax+1; grigliax+2].*(phi_y(ii)<=0);
% % %         % griglia = meshgrid(grigliax,grigliay);
% % %         
% % %         diagonale = 2*sqrt(2)*delta_x;
% % %         xPa = X(ki);
% % %         yPa = Y(ki);
% % %         modulo = sqrt(phi_x(ii).^2+phi_y(ii).^2);
% % %         xPb = X(ki)-diagonale.*phi_x(ii)./modulo;
% % %         yPb = Y(ki)-diagonale.*phi_y(ii)./modulo;
% % %         
% % %         for iter=1:(ceil(n^2/100))
% % %             xPh = (xPa+xPb)/2;
% % %             yPh = (yPa+yPb)/2;
% % %             Phia = interp2(X(griglia),Y(griglia),Phi(griglia),xPa,yPa);
% % %             Phib = interp2(X(griglia),Y(griglia),Phi(griglia),xPb,yPb);
% % %             Phih = interp2(X(griglia),Y(griglia),Phi(griglia),xPh,yPh);
% % %             if Phia*Phih<0
% % %                 xPb = xPh;
% % %                 yPb = yPh;
% % %             else
% % %                 xPa=xPh;
% % %                 yPa=yPh;
% % %             end
% % %         end
% % %         xC(ii)=xPh;
% % %         yC(ii)=yPh;
% % %     end
% % %     
% % %     % % % R = Phi(k);
% % %     % % % phi_x = (Phi(k+n)-Phi(k-n))/(2*delta_x);
% % %     % % % phi_y = (Phi(k+1)-Phi(k-1))/(2*delta_y);
% % % else
% % %     R = Phi(k);
% % %     modulo = sqrt(phi_x.^2+phi_y.^2);
% % %     xC = X(k) - R.*phi_x./modulo;
% % %     yC = Y(k) - R.*phi_y./modulo;
% % % end
% % % 
% % % x_C=xC;
% % % y_C=yC;
% % % xP = X(k);
% % % yP = Y(k);
% % % 
% % % quadrante = 1*(x_C>=xP & y_C>yP) + 2*(x_C<xP & y_C>yP) + 3*(x_C<xP & y_C<=yP) + 4*(x_C>=xP & y_C<=yP);
% % % 
% % % % R=
% % % % xC=x0+R*(xG-x0)/modulo;
% % % % yC=y0+R*(yG-y0)/modulo;