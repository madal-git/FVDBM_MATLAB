function [coord,result] = slice(cut,N_L,N_H,X1,X2,Y1,Y2,value_cell,value_nd,CELL,M)
% function [coord,result] = slice(cut,N_L,N_H,X1,X2,Y1,Y2,value_cell,value_nd)
% returns the series of point positions and the evaluated value at these
% points from a horizontal or vertical slice of the computation domain.
% cut is the geometric definition of the slice, has to be column vector
% N_L is the number of nodes on the horizontal outer surfaces
% N_H is the number of nodes on the vertical outer surfaces
% X1, X2, Y1 and Y2 is the boundaries of outer boundaries
% value_cell is the value from the centroid of triangles that is going to be
% sampled.
% value_nd is the value from the nodes that is going to be sampled.
% CELL is the cell data structure of all triangles
% M is the total number of triangles


%% reference size
e=((Y2-Y1)+(X2-X1))/2;
%% Check input data and creating coordinates of line of points on the slice
if cut(1,1)==Inf && cut(2,1)~=Inf
    msg='The slice is horizontal.';
    disp(msg);
    S=N_L;
    spc=(X2-X1)/(S-1);
    for i=1:S
        coord(1,i)=X1+(i-1)*spc;
        coord(2,i)=cut(2,1);
    end
elseif cut(1,1)~=Inf && cut(2,1)==Inf
    msg='The slice is vertical.';
    disp(msg);
    S=N_H;
    spc=(Y2-Y1)/(S-1);
    for i=1:S
        coord(2,i)=Y1+(i-1)*spc;
        coord(1,i)=cut(1,1);
    end
else
    error('The slice has to be vertical or horizontal!');
end

s1=length(value_cell(:,end));
s2=length(value_nd(:,end));
if s1~=s2
    error('The 2 sets of data that are going to be sampled are not the same type!');
end
l1=length(value_cell(1,:));
l2=length(value_nd(1,:));
if l1<l2
    error('The 2 sets of data that are going to be sampled are miss-located or the same type!');
end
result=zeros(s1,S);
%% Determin the value of two ends
% one end
[ND11,ND12,ratio1]=on_which_edge(coord(:,1),CELL,M);
if ND11==0 || ND12==0
    error('The intersection point on the outer boundary is not found');
else
    result(:,1)=(1-ratio1)*value_nd(:,ND11)+ratio1*value_nd(:,ND12);
end
% the other end
[ND21,ND22,ratio2]=on_which_edge(coord(:,S),CELL,M);
if ND21==0 || ND22==0
    error('The intersection point on the outer boundary is not found');
else
    result(:,S)=(1-ratio2)*value_nd(:,ND21)+ratio2*value_nd(:,ND22);
end
%% Filling the data in between
for i=1:S
    if i==1 || i==S
        ;
    else
        r=in_which_triangle(coord(:,i),CELL,M);
        if r==0 % The point is not within any triangle
            [nd1,nd2,ratio]=on_which_edge(coord(:,i),CELL,M);
            if nd1==0 || nd2==0 % The given point is not on any edge
                disp('The given point is out of computational domain!');
                result(:,i)=result(:,i-1);
            else
                result(:,i)=(1-ratio)*value_nd(:,nd1)+ratio*value_nd(:,nd2);
            end
        else % Found the triangle that enclose the point
            P=CELL{r};
            result(:,i)=value_cell(:,P{1}); % follow the assumption that the value in the triangle is constant
        end
    end
end