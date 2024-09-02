function d=dis(X1,X2)
% function d=dis(X1,X2) calculates the distance between two points
% X1 the first point, must be colum vector
% X2 the second point, must be column vector
% d is the distance between X1 an X2
% 
% N1=length(X1(1,:));
% M1=length(X1(:,1));
% 
% N2=length(X2(1,:));
% M2=length(X2(:,1));
% 
% if N1>1 || N2>1
%     error('The varibales must be a column vector!');
% end
% 
% if (M1<2 || M1>3) || (M2<2 || M2>3)
%     error('The given coornidates must be 2-D or 3-D');
% end
% 
% if M1~=M2
%     error('The given two points have different dimensions');
% end
% 
% M=M1;
% if M==2 % 2-D
    d=sqrt((X1(1,1)-X2(1,1))^2+(X1(2,1)-X2(2,1))^2);
% else % 3-D
%     d=sqrt((X1(1,1)-X2(1,1))^2+(X1(2,1)-X2(2,1))^2+(X1(3,1)-X2(3,1))^2);
% end