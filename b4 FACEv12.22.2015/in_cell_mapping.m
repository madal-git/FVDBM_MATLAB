function x=in_cell_mapping(FC,CELL,NODE,x_c,x_nd,s_cell,s_stcl_coord,s_map_id,s_map_coord,FMP,V)
% x=in_cell_mapping(FC,CELL,NODE,x_c,x_nd,s_cell,s_stcl_coord,s_map_id,s_map_coord,FMP,FX) returns 
% the mapped value at any stencil point.
% FC is the cell data structure of current face
% CELL is the enire cell data structure for CELL
% NODE is the entire cell data structure for NODE
% x_c is the entire any values data at centroids
% x_nd is the entire any value data at nodes
% s_cell is the vector storing the cell number that contains each stencil
% point
% s_stcl_coord is the coordinates of each stencil point
% s_map_id is the id number of three points that enclose each stencil point
% s_map_coord is the coordinates of the three points that enclose each stencil point
% FMP is the flag for different mapping method, FMP=1---First order;FMP=2---Second order

%%
if FMP==1
    x=x_c(:,s_cell);
elseif FMP==2
    if s_map_id(1,1)==0
        if s_map_id(2,1)==s_map_id(3,1) && s_map_id(2,1)==s_map_id(4,1) % The stencil point is located at a centroid
            if s_cell~=s_map_id(2,1)
                error('Wrong stencil!');
            end
            x=x_c(:,s_cell);
        elseif (s_map_id(2,1)~=s_map_id(3,1)) && (s_map_id(2,1)~=s_map_id(4,1)) && (s_map_id(3,1)~=s_map_id(4,1)) % a the stencil is located inside a triangle formed by three centroids
            if s_cell~=s_map_id(2,1) && (s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1))
                error('Wrong stencil!');
            end
            x=triangle_2ndmapping(C1,C2,C3,x_c(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
        else
            error('FC{20}, FC{28}, FC{35} and/or FC{42} contains false info!');
        end
    elseif s_map_id(1,1)==1
        if s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        C1=s_map_coord(1:2,1);
        C2=s_map_coord(3:4,1);
        C3=s_map_coord(5:6,1);
        x=triangle_2ndmapping(C1,C2,C3,x_nd(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
    elseif s_map_id(1,1)==2
        if s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        C1=s_map_coord(1:2,1);
        C2=s_map_coord(3:4,1);
        C3=s_map_coord(5:6,1);
        x=triangle_2ndmapping(C1,C2,C3,x_nd(:,s_map_id(2,1))',x_nd(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
    else
        error('FC{20} contains false info!');
    end
elseif FMP==3
    if s_map_id(1,1)==0
        if s_map_id(2,1)==s_map_id(3,1) && s_map_id(2,1)==s_map_id(4,1) % The stencil point is located at a centroid
            if s_cell~=s_map_id(2,1)
                error('Wrong stencil!');
            end
            x=x_c(:,s_cell);
        elseif (s_map_id(2,1)~=s_map_id(3,1)) && (s_map_id(2,1)~=s_map_id(4,1)) && (s_map_id(3,1)~=s_map_id(4,1)) % a the stencil is located inside a triangle formed by three centroids
            C1=s_map_coord(1:2,1);
            C2=s_map_coord(3:4,1);
            C3=s_map_coord(5:6,1);
            x=triangle_2ndmapping(C1,C2,C3,x_c(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
        else
            error('FC{20}, FC{28}, FC{35} and/or FC{42} contains false info!');
        end
    elseif s_map_id(1,1)==1
        if s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        C1=s_map_coord(1:2,1);
        C2=s_map_coord(3:4,1);
        C3=s_map_coord(5:6,1);
        x=triangle_2ndmapping(C1,C2,C3,x_nd(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
    elseif s_map_id(1,1)==2
        if s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        x=x_c(:,s_cell);
    else
        error('FC{20} contains false info!');
    end
elseif FMP==4
    if s_map_id(1,1)==0
        if s_map_id(2,1)==s_map_id(3,1) && s_map_id(2,1)==s_map_id(4,1) % The stencil point is located at a centroid
            if s_cell~=s_map_id(2,1)
                error('Wrong stencil!');
            end
            x=x_c(:,s_cell);
        elseif (s_map_id(2,1)~=s_map_id(3,1)) && ((s_map_id(2,1)~=s_map_id(4,1)) && (s_map_id(3,1)~=s_map_id(4,1))) % a the stencil is located inside a triangle formed by three centroids
            if s_cell~=s_map_id(2,1) && (s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1))
                error('Wrong stencil!');
            end
            C1=s_map_coord(1:2,1);
            C2=s_map_coord(3:4,1);
            C3=s_map_coord(5:6,1);
            x=triangle_2ndmapping(C1,C2,C3,x_c(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
        else
            error('FC{20}, FC{28}, FC{35} and/or FC{42} contains false info!');
        end
    elseif s_map_id(1,1)==1
        if s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        x=x_c(:,s_cell);
    elseif s_map_id(1,1)==2
        if s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        x=x_c(:,s_cell);
    else
        error('FC{20} contains false info!');
    end
elseif FMP==5
    if s_map_id(1,1)==0
        if s_map_id(2,1)==s_map_id(3,1) && s_map_id(2,1)==s_map_id(4,1) % The stencil point is located at a centroid
            if s_cell~=s_map_id(2,1)
                error('Wrong stencil!');
            end
            x=x_c(:,s_cell);
        elseif (s_map_id(2,1)~=s_map_id(3,1)) && (s_map_id(2,1)~=s_map_id(4,1)) && (s_map_id(3,1)~=s_map_id(4,1)) % a the stencil is located inside a triangle formed by three centroids
            x=x_c(:,s_cell);
        else
            error('FC{20}, FC{28}, FC{35} and/or FC{42} contains false info!');
        end
    elseif s_map_id(1,1)==1
        if s_cell~=s_map_id(3,1) && s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        C1=s_map_coord(1:2,1);
        C2=s_map_coord(3:4,1);
        C3=s_map_coord(5:6,1);
        x=triangle_2ndmapping(C1,C2,C3,x_nd(:,s_map_id(2,1))',x_c(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
    elseif s_map_id(1,1)==2
        if s_cell~=s_map_id(4,1)
            error('Wrong stencil!');
        end
        C1=s_map_coord(1:2,1);
        C2=s_map_coord(3:4,1);
        C3=s_map_coord(5:6,1);
        x=triangle_2ndmapping(C1,C2,C3,x_nd(:,s_map_id(2,1))',x_nd(:,s_map_id(3,1))',x_c(:,s_map_id(4,1))',s_stcl_coord)';
    else
        error('FC{20} contains false info!');
    end
else
    error('The flag for mapping method is incorrect!');
end