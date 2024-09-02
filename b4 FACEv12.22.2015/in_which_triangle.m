function x=in_which_triangle(p,CELL,M)
% function x=in_which_triangle(p,CELL,M) returns the triangle number that
% enclose the given point. if none of the triangles containes the given
% point, then return zero.
% p is the coordinate of given point, has to be column vector
% CELL is the cell data structure of all triangles
% M is the total number of triangles

for r=1:M
    P=CELL{r};
    if in_triangle(p,P{13},P{14},P{15})
        break;
    end
end

% check
if r==M
    if in_triangle(p,P{13},P{14},P{15})
        x=r;
    else
        x=0;
    end
else
    x=r;
end