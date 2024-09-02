function b = on_which_boundary(C,X1,X2,Y1,Y2)
% function b = on_which_boundary(C,X1,X2,Y1,Y2) determines on which bounday
% of of rectangular domain the given point is located
% b is the return value for the boundary indicator
%    b=0---not on any boundary
%    b=1---on top boundary
%    b=2---on right boundary
%    b=3---on bottom boundary
%    b=4---on left boundary
% C is the coordinate of the giben point, a column vector
% X1, X2, Y1 and Y2 are the bounds of the rectangular domain

BL=[X1;Y1];
BR=[X2;Y1];
TL=[X1;Y2];
TR=[X2;Y2];

on_top=on_edge(C,TL,TR);
on_right=on_edge(C,TR,BR);
on_bottom=on_edge(C,BL,BR);
on_left=on_edge(C,BL,TL);

if (on_top+on_right+on_bottom+on_left)==0 % Not on any outer boundary
    b=0;
elseif (on_top+on_right+on_bottom+on_left)==1 % On one of the outer boundary
    if on_top
        b=1;
    elseif on_right
        b=2;
    elseif on_bottom
        b=3;
    elseif on_left
        b=4;
    else
        error('Logical error!');
    end
elseif (on_top+on_right+on_bottom+on_left)==2 % At one of the corner points shared by two boundaries
    if on_top && on_right
        b=[1,2];
    elseif on_right && on_bottom
        b=[2,3];
    elseif on_bottom && on_left
        b=[3,4];
    elseif on_left && on_top
        b=[4,1];
    else
        error('Logical error!');
    end
else
    error('Logical error!');
end