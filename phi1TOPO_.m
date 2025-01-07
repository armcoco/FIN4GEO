function RESULT = phi1TOPO_(x,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traccia=load('profilo50.dat');

[l,m]=size(x);
 
for i=1:l
    for j=1:m
        
        if abs(x(i,j))<40000
            dist=(x(i,j)-traccia(:,1)).^2+(y(i,j)-traccia(:,2)).^2;
            [val,ind]=sort(dist);
            a=(traccia(ind(2),2)-traccia(ind(1),2))/(traccia(ind(2),1)-traccia(ind(1),1));
            b=-1;
            c=traccia(ind(1),2)-a*traccia(ind(1),1);   
            RESULT(i,j)=abs(a*x(i,j)+b*y(i,j)+c)./(sqrt(a.^2+b.^2));
            
            if ind(1)<ind(2)
                t=[traccia(ind(2),1)-traccia(ind(1),1);traccia(ind(2),2)-traccia(ind(1),2)];
                d=[x(i,j)-traccia(ind(1),1);y(i,j)-traccia(ind(1),2)];	
            else
                t=[traccia(ind(1),1)-traccia(ind(2),1);traccia(ind(1),2)-traccia(ind(2),2)];
                d=[x(i,j)-traccia(ind(2),1);y(i,j)-traccia(ind(2),2)];
            end
            p=t(1)*d(2)-d(1)*t(2);
            if p<0
                RESULT(i,j)=-RESULT(i,j);
            end
            
      elseif x(i,j)<=-40000
            RESULT(i,j)=y(i,j)-traccia(1,2);
         else
             RESULT(i,j)=y(i,j)-traccia(end,2);    
        end
        
    end
end

RESULT(:,1)=RESULT(:,2);
RESULT(:,end)=RESULT(:,end-1);
RESULT(1,:)=RESULT(2,:);
RESULT(end,:)=RESULT(end-1,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Profilo Reale Etna %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% traccia0=load('profilo50.dat');
% traccia=interp1(traccia0(:,1),traccia0(:,2),x(1,:));
% i=find(isnan(traccia));
% traccia(i)=0;
% traccia=[x(1,:)' traccia'];
%  
% [l,m]=size(x);
% D=zeros(l,m);
% 
% for i=1:l
%     for j=1:m
%         
%         if abs(x(i,j))<40000
%             dist=abs(x(i,j)-traccia(:,1));
%             [val,ind]=sort(dist);
%             RESULT(i,j)=y(i,j)-traccia(ind(1),2);           
%         elseif x(i,j)<=-40000
%             RESULT(i,j)=y(i,j)-traccia0(1,2);
%         else
%             RESULT(i,j)=y(i,j)-traccia0(end,2);    
%         end
%         
%     end
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% END Profilo Reale Etna %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% %for test
% 
% global c_parameter m_par;
% c_parameter=10000;
% m_par=0.5
% tank=[-1 1 -1 1]; %*40*scale;
% 
% a=tank(1);
% b=tank(2);
% c=tank(3);
% d=tank(4);
% 
% m=300;
% n=300;
% dx=(b-a)/m; %spatial step in the x-diricetion
% xi=a:dx:b;
% dy=(d-c)/n; %spatial step in the y-diricetion
% yi=c:dy:d;
% 
% [X,Y] = meshgrid(xi,yi);
% [nn,mm]=size(X);
% 
% x=xt_(X);
% y=yt_(Y); 
% 
% traccia0=load('profilo50.dat');
% traccia=interp1(traccia0(:,1),traccia0(:,2),x(1,:));
% i=find(isnan(traccia));
% traccia(i)=0;
% traccia=[x(1,:)' traccia'];
%  
% [l,m]=size(x);
% D=zeros(l,m);
% 
% for i=1:l
%     for j=1:m
%         
%         if abs(x(i,j))<40000
%             dist=(x(i,j)-traccia(:,1)).^2+(y(i,j)-traccia(:,2)).^2;
%             [val,ind]=sort(dist);
%             RESULT(i,j)=y(i,j)-traccia(ind(1),2);           
%         else
%             RESULT(i,j)=y(i,j);
%         end
%         
%     end
% end
% 
% RESULT(:,1)=RESULT(:,2);
% RESULT(:,end)=RESULT(:,end-1);
% 
% 
% 
% %OLD Routine funzionante quando griglia uniforme 
% %  
% % for i=1:l
% %     for j=1:m
% %         
% %         if abs(x(i,j))<40000
% %             dist=(x(i,j)-traccia(:,1)).^2+(y(i,j)-traccia(:,2)).^2;
% %             [val,ind]=sort(dist);
% %             a=(traccia(ind(2),2)-traccia(ind(1),2))/(traccia(ind(2),1)-traccia(ind(1),1));
% %             b=-1;
% %             c=traccia(ind(1),2)-a*traccia(ind(1),1);   
% %             RESULT(i,j)=abs(a*x(i,j)+b*y(i,j)+c)./(sqrt(a.^2+b.^2));
% %             
% %             if ind(1)<ind(2)
% %                 t=[traccia(ind(2),1)-traccia(ind(1),1);traccia(ind(2),2)-traccia(ind(1),2)];
% %                 d=[x(i,j)-traccia(ind(1),1);y(i,j)-traccia(ind(1),2)];	
% %             else
% %                 t=[traccia(ind(1),1)-traccia(ind(2),1);traccia(ind(1),2)-traccia(ind(2),2)];
% %                 d=[x(i,j)-traccia(ind(2),1);y(i,j)-traccia(ind(2),2)];
% %             end
% %             p=t(1)*d(2)-d(1)*t(2);
% %             if p<0
% %                 RESULT(i,j)=-RESULT(i,j);
% %             end
% %             
% %         else
% %             RESULT(i,j)=y(i,j);
% %         end
% %         
% %     end
% % end
% % 
% % RESULT(:,1)=RESULT(:,2);
% % RESULT(:,end)=RESULT(:,end-1);
% % 
% 
% % for test
% [c,h]=contour(x,y,RESULT);
% clabel(c,h)
% hold on
% plot(traccia(:,1),traccia(:,2),'o')
% 
% aaa=reshape(x,prod(size(x)),1);
% size(aaa)
% scatter(reshape(x,prod(size(x)),1),reshape(y,prod(size(x)),1),1,reshape(RESULT,prod(size(x)),1))
% hold off
% 
% figure
% plot(RESULT)
% [a,b]=find(isnan(RESULT))
% x(end,end)
% y(end,end)