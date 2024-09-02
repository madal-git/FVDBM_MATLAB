function y = inter_extrp(x1,x2,x3,y1,y2,y3,x,fie)

% fucntion y = inter_extrp(x1,x2,x3,y1,y2,y3,x,fie) linearly or cubically
% interpolate or extrapolate the value y at location x with known locations
% x1, x2, and x3 and their corresponding values y1, y2, and y3
% fie is the controller for linear or cubic scheme 
% fie=0----linear interpolation and extrapolation, x3 and y3 will be
% ignored
% fie=1----cubic interpolation and extrapolation
[r1,c1]=size(y1);
[r2,c2]=size(y2);
[r3,c3]=size(y3);
if (r1~=r2 || r2~=r3) || (c1~=c2 || c2~=c3)
    error('the input values must have the same dimensions!');
end
y=zeros(r1,c1);
if fie==0
    A=[x1, 1; x2, 1];
    for i=1:r1
        for j=1:c1
            B=[y1(i,j); y2(i,j)];
            y(i,j)=[x, 1]*(A\B);
        end
    end
elseif fie==1
    A=[x1^2, x1, 1; x2^2, x2, 1; x3^2, x3, 1];
    for i=1:r1
        for j=1:c1
            B=[y1(i,j); y2(i,j); y3(i,j)];
            y(i,j)=[x^2, x, 1]*(A\B);
        end
    end
else
    error('Wrong flag for interpolation/extrapolation!');
end