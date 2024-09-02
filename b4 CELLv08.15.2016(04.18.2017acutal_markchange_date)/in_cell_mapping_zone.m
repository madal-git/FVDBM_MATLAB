function x=in_cell_mapping_zone(P,Zone,x_centroid,x_nd1,x_nd2,x_nd3,C,FMP)
% function x=in_cell_mapping_b(P,Zone,x_centroid,x_nd1,x_nd2,x_nd3,C,FMP) returns 
% the mapped value at any point located inside of a cell that is attached to boundary.
% P is the cell data of current cell, P=CELL{r}
% Zone is the Zone ID of the target point that is located in the current
   % cell P
% x_centroid is the value at the centroid of the cell
% x_nd1 is the value at the first node of the cell
% x_nd2 is the value at the second node of the cell
% x_nd3 is the value at the third node of the cell
% C is the coordinate of the target point, column vector
% FMP is the flag for different mapping method, FMP=1---First order;FMP=2---Second order


if FMP==1
    x=x_centroid;
elseif FMP==2
    if Zone==4
        x=x_centroid;
    elseif Zone==1
        C1=P{5};
        C2=P{13};
        C3=P{14};
        x=triangle_2ndmapping(C1,C2,C3,x_centroid',x_nd1',x_nd2',C)';
    elseif Zone==2
        C1=P{5};
        C2=P{14};
        C3=P{15};
        x=triangle_2ndmapping(C1,C2,C3,x_centroid',x_nd2',x_nd3',C)';
    elseif Zone==3
        C1=P{5};
        C2=P{15};
        C3=P{13};
        x=triangle_2ndmapping(C1,C2,C3,x_centroid',x_nd3',x_nd1',C)';
    else
        error('The zone info is incorrect!');
    end
%     x=x_centroid;
else
    error('The flag for mapping method is incorrect!');
end