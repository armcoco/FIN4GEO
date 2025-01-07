
MS=20;
obj=2;
while 1==1
    [xg,yg,button]=ginput(1);
    if any(ishandle(obj))
        delete(obj);
    end
    j=ceil((xg-a)/dx);
    i=ceil((yg-c)/dy);
    i=max(i,1);
    j=max(j,1);
    K_temp=i+(j-1)*(nn);
    KK=[K_temp, K_temp+1, K_temp+nn, K_temp+1+nn];
    [~,Ki]=min((X(KK)-xg).^2+(Y(KK)-yg).^2);
    K=KK(Ki);
    stencil_temp=find(abs(M(K,:))>1e-10);
    st1=stencil_temp(stencil_temp<=nn*mm);
    st2=stencil_temp(stencil_temp>nn*mm)-nn*mm;
    if button==1
        obj=plot(X(st1),Y(st1),'ob','markerSize',MS);
        obj(2)=plot(X(st2),Y(st2),'dr','markerSize',MS);
    elseif button==3
        break;
    end
end
