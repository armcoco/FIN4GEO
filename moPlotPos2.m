function AxisPos = moPlotPos2(nCols, nRows, Rect,factorW,factorH)

if nargin<=3
    factorW=0.03;
    factorH=0.03;
end
AxisPos=zeros(nCols*nRows,4);
Width_off=factorW*Rect(3)/nCols;
Width=(Rect(3)-(nCols-1)*Width_off)/nCols;
Height_off=factorH*Rect(4)/nRows;
Height=(Rect(4)-(nRows-1)*Height_off)/nRows;
AxisPos(:,3)=Width;
AxisPos(:,4)=Height;
X_temp=repmat([Rect(1)+(0:(nCols-1))*(Width+Width_off)],nRows,1);
AxisPos(:,1)=X_temp(:);
AxisPos(:,2)=repmat([Rect(2)+(0:(nRows-1))*(Height+Height_off)],1,nCols);
