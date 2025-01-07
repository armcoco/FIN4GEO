function RESULT = phi1Topo(x,y,filename)

% if topo==0
%     RESULT = y; %topo flat
% elseif topo==1
%     a=3000;
%     b=3000;
%     RESULT= y-a^2*b./(x.^2+a^2);%topo gaussiana: a^2*b/(r^2+a^2)
% else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Profilo Gaussiano%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a=1600;
% b=0;
% c=6000;
% RESULT=y-(a*exp(-(x-b).^2/c^2)-1300);
% return 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% Profilo Reale Vulcano %%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return
%%%for test
% [x,y]=meshgrid(0:100:20000,-4000:100:5000); 

traccia0=load(filename);
 
[l,m]=size(x);
X=(x(1,1:end-1)+x(1,2:end))/2;
traccia=interp1(traccia0(:,1),traccia0(:,2),X');
traccia=[X' traccia];
aa=find(traccia(:,1)>15000);
traccia(aa,2)=traccia0(1,2);

% figure(100)
% plot(traccia(:,1),traccia(:,2))
 
for i=1:l
    for j=1:m
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
        
    end
end
RESULT(:,end)=RESULT(:,end-1);
RESULT(end,:)=RESULT(end-1,:);
RESULT(1,:)=RESULT(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END Profilo Reale Etna %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end

% % % for test
% [c,h]=contour(x,y,RESULT);
% clabel(c,h)
% hold on
% plot(traccia(:,1),traccia(:,2),'o')
% hold off
