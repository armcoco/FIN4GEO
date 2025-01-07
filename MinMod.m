function MM=MinMod(A,B)

%n=length(A(:));

if (size(A)~=size(B))
    disp('ERROR: size of A and size of B must be equal');
end

MM=zeros(size(A));
mask=find(A.*B>0);

MM(mask)=A(mask).*(abs(A(mask))<=abs(B(mask)))+B(mask).*(abs(A(mask))>abs(B(mask)));

% for i=1:n
%     a=A(i);
%     b=B(i);
%     
%     
%     if (abs(a)<=abs(b) && a*b>0)
%         mm(i)=a;
%     elseif (abs(a)>abs(b) && a*b>0)
%         mm(i)=b;
%     elseif (a*b<=0)
%         mm(i)=0;
%     end
% end
% mm=mm';
