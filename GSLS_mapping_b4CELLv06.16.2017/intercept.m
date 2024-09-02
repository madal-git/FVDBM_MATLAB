function [Bool_found, C_i] = intercept(Ori,n,C1,C2)
%% function [Bool_found, C_i] = intercept(Ori,n,C1,C2) finds the coordinates of intercept point
% between a segment bounded by C1, C2 and a vector starting from a point named Ori with a direction n

% Bool_found is the boolean variable indicating whether the intercept point
% is found or not.
% C_i is the coordinates of the intercept point, a column vector
% Ori is the coordinates of the origin point of the vector, a column vector
% n is the direction of the vector, a column vector, which does not have to
% be a unit vector
% C1 is the coordinates of the first end point of the segment, a column
% vector
% C2 is the coordinates of the second end point of the segment, a column
% vector
Tol0=1e-10;
Tol1=1e-6;
Tol2=1e-7;
Tol3=1e-9; % This value is from 1e-10 to 3e-10

% Check
if dis(C1,C2)<Tol2
    error('The segment to be intercepted is too short!');
end
if norm(n)<Tol3
    error('The vector is a zero vector!');
end
v1=C1-Ori;
v2=C2-Ori;
v_mid=(C1+C2)/2-Ori;
if n'*v_mid<0
    n=-n; % This is due to asin() could only return a angle between 0 and 90 degrees
end
ang1=angle(v1,n);
ang2=angle(v2,n);
ang1_2=angle(v1,v2);

%% Algorithm 1
% if single(Tol0+ang1+ang2)==single(Tol0+ang1_2)
%     Bool_found=1;
%     if single(Tol1+ang1)==single(Tol1)
%         C_i=C1;
%     elseif single(Tol1+ang2)==single(Tol1)
%         C_i=C2;
%     else
%         if abs(C1(1,1)-C2(1,1))<Tol2 % The segment to be intercepted is vertical
%             x=C1(1,1);
%             if abs(n(1,1)/norm(n))<Tol2
%                 error('Logic error!');
%             end
%             ratio=(x-Ori(1,1))/n(1,1);
%             y=Ori(2,1)+ratio*n(2,1);
%         elseif abs(C1(2,1)-C2(2,1))<Tol2 % The segment to be intercepted is horizantal
%             y=C1(2,1);
%             if abs(n(2,1)/norm(n))<Tol2
%                 error('Logic error!');
%             end
%             ratio=(y-Ori(2,1))/n(2,1);
%             x=Ori(1,1)+ratio*n(1,1);
%         else % The segment to be intercepted is a regular slope
%             % Construct y=a1*x+b1 from C1 and C2
%             if C1(1,1)<C2(1,1)
%                 a1=(C2(2,1)-C1(2,1))/(C2(1,1)-C1(1,1));
%             else
%                 a1=(C1(2,1)-C2(2,1))/(C1(1,1)-C2(1,1));
%             end
%             b1_1=C1(2,1)-a1*C1(1,1);
%             b1_2=C2(2,1)-a1*C2(1,1);
%             if abs(b1_1-b1_2)>Tol3
%                 error('Logic error!');
%             end
%             b1=b1_1;
%             if abs(n(1,1)/norm(n))<Tol2 % The vector that is going to intercept the segment is vertical
%                 x=Ori(1,1);
%                 y=a1*x+b1;
%             elseif abs(n(2,1)/norm(n))<Tol2 % The vector that is going to intercept the segment is horizontal
%                 y=Ori(2,1);
%                 x=(y-b1)/a1;
%             else % The vector that is going to intercept the segment is a regular slope
%                 % Construct y=a2*x+b2 from n and Ori
%                 a2=n(2,1)/n(1,1);
%                 b2=Ori(2,1)-a2*Ori(1,1);
%                 if abs(a1-a2)<Tol2
%                     error('Logic error!');
%                 end
%                 x=(b2-b1)/(a1-a2);
%                 y=(a1*b2-a2*b1)/(a1-a2);
%             end
%         end
%         C_i=[x;y];
%     end
% else
%     Bool_found=0;
%     C_i=[inf;inf];
% end

%% Algorithm 2
ang=angle(n,(C2-C1));
if ang<Tol1 || abs(ang-180)<Tol1 % Parallel
    Bool_found=0;
    C_i=[inf;inf];
else
    if abs(C1(1,1)-C2(1,1))/dis(C1,C2)<Tol2 % The segment to be intercepted is vertical
        x=(C1(1,1)+C2(1,1))/2;
        if abs(n(1,1)/norm(n))<Tol2
            error('Logic error!');
        end
        ratio=(x-Ori(1,1))/n(1,1);
        y=Ori(2,1)+ratio*n(2,1);
%         y=Ori(2,1);
    elseif abs(C1(2,1)-C2(2,1))/dis(C1,C2)<Tol2 % The segment to be intercepted is horizantal
        y=(C1(2,1)+C2(2,1))/2;
        if abs(n(2,1)/norm(n))<Tol2
            error('Logic error!');
        end
        ratio=(y-Ori(2,1))/n(2,1);
        x=Ori(1,1)+ratio*n(1,1);
%         x=Ori(1,1);
    else % The segment to be intercepted is a regular slope
        % Construct y=a1*x+b1 from C1 and C2
        if C1(1,1)<C2(1,1)
            a1=(C2(2,1)-C1(2,1))/(C2(1,1)-C1(1,1));
        else
            a1=(C1(2,1)-C2(2,1))/(C1(1,1)-C2(1,1));
        end
        b1_1=C1(2,1)-a1*C1(1,1);
        b1_2=C2(2,1)-a1*C2(1,1);
        if abs(b1_1-b1_2)/dis(C1,C2)>Tol3
            error('Logic error!');
        end
        b1=(b1_1+b1_2)/2;
        if abs(n(1,1)/norm(n))<Tol2 % The vector that is going to intercept the segment is vertical
            x=Ori(1,1);
            y=a1*x+b1;
        elseif abs(n(2,1)/norm(n))<Tol2 % The vector that is going to intercept the segment is horizontal
            y=Ori(2,1);
            x=(y-b1)/a1;
        else % The vector that is going to intercept the segment is a regular slope
            % Construct y=a2*x+b2 from n and Ori
            a2=n(2,1)/n(1,1);
            b2=Ori(2,1)-a2*Ori(1,1);
            if abs(a1-a2)<Tol2
                error('Logic error!');
            end
            x=(b2-b1)/(a1-a2);
            y=(a1*b2-a2*b1)/(a1-a2);
        end
    end
    C_i=[x;y];
    if abs(dis(C_i,C1)+dis(C_i,C2)-dis(C1,C2))/dis(C1,C2)<Tol2
        Bool_found=1;
        if dis(C_i,C1)<Tol2
            C_i=C1;
        elseif dis(C_i,C2)<Tol2
            C_i=C2;
        else
            ;
        end
    else
        Bool_found=0;
        C_i=[inf;inf];
    end
end
