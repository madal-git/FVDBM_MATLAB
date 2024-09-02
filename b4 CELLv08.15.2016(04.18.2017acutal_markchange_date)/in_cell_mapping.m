function x=in_cell_mapping(FC,CELL,NODE,x_c,x_nd,s_cell,s_stcl_coord,s_map_id,s_map_coord,gradient,FMP,V)
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
% gradient is the gradient of the PDFs at each cell centroid
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
elseif FMP==6 % Least squre scheme
    %% Algorithm 1---no upwind
    cell_center=CELL{s_cell};
    n=s_stcl_coord-cell_center{5};
    x=x_c(:,s_cell)+gradient(:,:,s_cell)'*n;
    %% Algorithm 2---upwind is determined by V and s_stcl_coord-cell_center{5}, and the upwind components take the ls value(2nd order), the downwind ones take the centroid value(1st order)
%     cell_center=CELL{s_cell};
%     n=s_stcl_coord-cell_center{5};
%     x_pool=[x_c(:,s_cell),x_c(:,s_cell)+gradient(:,:,s_cell)'*n];
%     upwind_pointer=V'*n;
%     upwind_pointer=(upwind_pointer>=0);
%     upwind_pointer_pool=[~upwind_pointer,upwind_pointer];
%     x_pool_upwind=x_pool.*upwind_pointer_pool;
%     x=x_pool_upwind(:,1)+x_pool_upwind(:,2);
    %% Algorithm 3---upwind is determined by V and s_stcl_coord-cell_center{5}, and the upwind components take the centroid value(1st order), the downwind ones take the ls value(2nd order)
%     cell_center=CELL{s_cell};
%     n=s_stcl_coord-cell_center{5};
%     x_pool=[x_c(:,s_cell),x_c(:,s_cell)+gradient(:,:,s_cell)'*n];
%     upwind_pointer=V'*n;
%     upwind_pointer=(upwind_pointer>=0);
%     upwind_pointer_pool=[upwind_pointer,~upwind_pointer];
%     x_pool_upwind=x_pool.*upwind_pointer_pool;
%     x=x_pool_upwind(:,1)+x_pool_upwind(:,2);
    %% Algorithm 4---upwind is determined by V and the direction of stencil, and the upwind components take the ls value(2nd order), the downwind ones take the centroid value(1st order)
%     cell_center=CELL{s_cell};
%     n=FC{4};
%     x_pool=[x_c(:,s_cell),x_c(:,s_cell)+gradient(:,:,s_cell)'*(s_stcl_coord-cell_center{5})];
%     upwind_pointer=(n*V)';
%     upwind_pointer=(upwind_pointer>0);
%     upwind_pointer_pool=[~upwind_pointer,upwind_pointer];
%     x_pool_upwind=x_pool.*upwind_pointer_pool;
%     x=x_pool_upwind(:,1)+x_pool_upwind(:,2);
    %% Algorithm 5---upwind is determined by V and the direction of stencil, and the upwind components take the centroid value(1st order), the downwind ones take the ls value(2nd order)
%     cell_center=CELL{s_cell};
%     n=FC{4};
%     x_pool=[x_c(:,s_cell),x_c(:,s_cell)+gradient(:,:,s_cell)'*(s_stcl_coord-cell_center{5})];
%     upwind_pointer=(n*V)';
%     upwind_pointer=(upwind_pointer>0);
%     upwind_pointer_pool=[upwind_pointer,~upwind_pointer];
%     x_pool_upwind=x_pool.*upwind_pointer_pool;
%     x=x_pool_upwind(:,1)+x_pool_upwind(:,2);
else
    error('The flag for mapping method is incorrect!');
end