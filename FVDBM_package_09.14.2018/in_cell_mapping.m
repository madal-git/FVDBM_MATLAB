function x=in_cell_mapping(FC,CELL,NODE,x_c,x_nd,s_cell,s_stcl_coord,s_map_id,s_map_coord,gradient,FMP,V,FPDC)
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
% FPDC is the flag for periodic boundary conditions. FPDC=0---No periodic
% boundaries; FPDC=1---Only left & right boundaries are periodic;
% FPDC=2---Only top & bottom boundaries are periodic; FPDC=3---All
% boundaries are periodic

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
elseif FMP==7 % Xiaofeng's scheme
    cell_center=CELL{s_cell};
    %% Using three point, one of which is from boundary
    if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
        cell_neighbor=cell_center{38};
        CoeA=cell_center{40};
        F=[(x_c(:,cell_neighbor(1))-x_c(:,s_cell))';(x_c(:,cell_neighbor(2))-x_c(:,s_cell))';(x_c(:,cell_neighbor(3))-x_c(:,s_cell))'];
%         x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
        x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
    elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
        cell_neighbor=cell_center{41};
        fc_nd=cell_center{42};
        CoeA=cell_center{44};
        if cell_neighbor(1)==0
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_1=2*fs_intcpt-x_c(:,s_cell);
%             x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
            x_c_2=x_c(:,cell_neighbor(2));
            x_c_3=x_c(:,cell_neighbor(3));
        elseif cell_neighbor(2)==0
            x_c_1=x_c(:,cell_neighbor(1));
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_2=2*fs_intcpt-x_c(:,s_cell);
%             x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
            x_c_3=x_c(:,cell_neighbor(3));
        elseif cell_neighbor(3)==0
            x_c_1=x_c(:,cell_neighbor(1));
            x_c_2=x_c(:,cell_neighbor(2));
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_3=2*fs_intcpt-x_c(:,s_cell);
%             x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
        else
            error('logic error!');
        end
        F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%         x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
        x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
%         x=x_c(:,s_cell);
    elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
        if FPDC==0 % No periodic boundaries
            cell_neighbor=cell_center{41};
            fc_nd=cell_center{42};
            CoeA=cell_center{44};
            if cell_neighbor(1)==0
%                 x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_1=2*fs_intcpt-x_c(:,s_cell);
                x_c_2=x_c(:,cell_neighbor(2));
                x_c_3=x_c(:,cell_neighbor(3));
            elseif cell_neighbor(2)==0
                x_c_1=x_c(:,cell_neighbor(1));
%                 x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_2=2*fs_intcpt-x_c(:,s_cell);
                x_c_3=x_c(:,cell_neighbor(3));
            elseif cell_neighbor(3)==0
                x_c_1=x_c(:,cell_neighbor(1));
                x_c_2=x_c(:,cell_neighbor(2));
%                 x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_3=2*fs_intcpt-x_c(:,s_cell);
            else
                error('logic error!');
            end
            F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
            x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
        elseif FPDC==1 % Only left & right boundaries are periodic
            cell_neighbor=cell_center{49};
            fc_nd=cell_center{50};
            CoeA=cell_center{52};
            len=length(setxor(cell_neighbor,0));
            if len==2 % The current cell is attached to a non-periodic boundary
                if cell_neighbor(1)==0
%                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_1=2*fs_intcpt-x_c(:,s_cell);
                    x_c_2=x_c(:,cell_neighbor(2));
                    x_c_3=x_c(:,cell_neighbor(3));
                elseif cell_neighbor(2)==0
                    x_c_1=x_c(:,cell_neighbor(1));
%                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_2=2*fs_intcpt-x_c(:,s_cell);
                    x_c_3=x_c(:,cell_neighbor(3));
                elseif cell_neighbor(3)==0
                    x_c_1=x_c(:,cell_neighbor(1));
                    x_c_2=x_c(:,cell_neighbor(2));
%                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_3=2*fs_intcpt-x_c(:,s_cell);
                else
                    error('logic error!');
                end
            elseif len==4 % The current cell is attached to a periodic boundary
                x_c_1=x_c(:,cell_neighbor(1));
                x_c_2=x_c(:,cell_neighbor(2));
                x_c_3=x_c(:,cell_neighbor(3));
            else
                error('logic error!');
            end
            F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
            x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
        elseif FPDC==2 % Only top & bottom boundaries are periodic
            cell_neighbor=cell_center{53};
            fc_nd=cell_center{54};
            CoeA=cell_center{56};
            len=length(setxor(cell_neighbor,0));
            if len==2 % The current cell is attached to a non-periodic boundary
                if cell_neighbor(1)==0
%                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_1=2*fs_intcpt-x_c(:,s_cell);
                    x_c_2=x_c(:,cell_neighbor(2));
                    x_c_3=x_c(:,cell_neighbor(3));
                elseif cell_neighbor(2)==0
                    x_c_1=x_c(:,cell_neighbor(1));
%                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_2=2*fs_intcpt-x_c(:,s_cell);
                    x_c_3=x_c(:,cell_neighbor(3));
                elseif cell_neighbor(3)==0
                    x_c_1=x_c(:,cell_neighbor(1));
                    x_c_2=x_c(:,cell_neighbor(2));
%                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                    fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                    x_c_3=2*fs_intcpt-x_c(:,s_cell);
                else
                    error('logic error!');
                end
%                 x=x_c(:,s_cell);
            elseif len==4 % The current cell is attached to a periodic boundary
                x_c_1=x_c(:,cell_neighbor(1));
                x_c_2=x_c(:,cell_neighbor(2));
                x_c_3=x_c(:,cell_neighbor(3));
            else
                error('logic error!');
            end
            F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
            x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';   
        elseif FPDC==3 % All boundaries are periodic
            cell_neighbor=cell_center{45};
            CoeA=cell_center{48};
            F=[(x_c(:,cell_neighbor(1))-x_c(:,s_cell))';(x_c(:,cell_neighbor(2))-x_c(:,s_cell))';(x_c(:,cell_neighbor(3))-x_c(:,s_cell))'];
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
            x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
        else
            error('Logic error!');
        end
    else
        error('Logic error!');
    end
    %% Using only existing centroids, don't use boundary nodes at all
%     if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{38};
%         CoeA=cell_center{40};
%         F=[(x_c(:,cell_neighbor(1))-x_c(:,s_cell))';(x_c(:,cell_neighbor(2))-x_c(:,s_cell))';(x_c(:,cell_neighbor(3))-x_c(:,s_cell))'];
% %         x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%         x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
%     elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{41};
%         fc_nd=cell_center{42};
%         CoeA=cell_center{44};
%         if cell_neighbor(1)==0
%             x_c_2=x_c(:,cell_neighbor(2));
%             x_c_3=x_c(:,cell_neighbor(3));
%             F=[(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%         elseif cell_neighbor(2)==0
%             x_c_1=x_c(:,cell_neighbor(1));
%             x_c_3=x_c(:,cell_neighbor(3));
%             F=[(x_c_1-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%         elseif cell_neighbor(3)==0
%             x_c_1=x_c(:,cell_neighbor(1));
%             x_c_2=x_c(:,cell_neighbor(2));
%             F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))'];
%         else
%             error('logic error!');
%         end
%         
%         x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
%     elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
%         if FPDC==0 % No periodic boundaries
%             cell_neighbor=cell_center{41};
%             fc_nd=cell_center{42};
%             CoeA=cell_center{44};
%             if cell_neighbor(1)==0
% %                 x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                 x_c_1=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                 x_c_2=x_c(:,cell_neighbor(2));
%                 x_c_3=x_c(:,cell_neighbor(3));
%                 F=[(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             elseif cell_neighbor(2)==0
%                 x_c_1=x_c(:,cell_neighbor(1));
% %                 x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                 x_c_2=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                 x_c_3=x_c(:,cell_neighbor(3));
%                 F=[(x_c_1-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             elseif cell_neighbor(3)==0
%                 x_c_1=x_c(:,cell_neighbor(1));
%                 x_c_2=x_c(:,cell_neighbor(2));
% %                 x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                 x_c_3=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                 F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))'];
%             else
%                 error('logic error!');
%             end
% %             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
% %             x=x_c(:,s_cell);
%         elseif FPDC==1 % Only left & right boundaries are periodic
%             cell_neighbor=cell_center{49};
%             fc_nd=cell_center{50};
%             CoeA=cell_center{52};
%             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
%                 if cell_neighbor(1)==0
% %                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_1=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     x_c_2=x_c(:,cell_neighbor(2));
%                     x_c_3=x_c(:,cell_neighbor(3));
%                     F=[(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%                 elseif cell_neighbor(2)==0
%                     x_c_1=x_c(:,cell_neighbor(1));
% %                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_2=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     x_c_3=x_c(:,cell_neighbor(3));
%                     F=[(x_c_1-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%                 elseif cell_neighbor(3)==0
%                     x_c_1=x_c(:,cell_neighbor(1));
%                     x_c_2=x_c(:,cell_neighbor(2));
% %                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_3=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))'];
%                 else
%                     error('logic error!');
%                 end
% %                 x=x_c(:,s_cell);
%             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
%                 x_c_1=x_c(:,cell_neighbor(1));
%                 x_c_2=x_c(:,cell_neighbor(2));
%                 x_c_3=x_c(:,cell_neighbor(3));
%                 F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             else
%                 error('logic error!');
%             end
%             
% %             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
%         elseif FPDC==2 % Only top & bottom boundaries are periodic
%             cell_neighbor=cell_center{53};
%             fc_nd=cell_center{54};
%             CoeA=cell_center{56};
%             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
%                 if cell_neighbor(1)==0
% %                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_1=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     x_c_2=x_c(:,cell_neighbor(2));
%                     x_c_3=x_c(:,cell_neighbor(3));
%                     F=[(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%                 elseif cell_neighbor(2)==0
%                     x_c_1=x_c(:,cell_neighbor(1));
% %                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_2=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     x_c_3=x_c(:,cell_neighbor(3));
%                     F=[(x_c_1-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%                 elseif cell_neighbor(3)==0
%                     x_c_1=x_c(:,cell_neighbor(1));
%                     x_c_2=x_c(:,cell_neighbor(2));
% %                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
%                     x_c_3=x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2))-x_c(:,s_cell);
%                     F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))'];
%                 else
%                     error('logic error!');
%                 end
% %                 x=x_c(:,s_cell);
%             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
%                 x_c_1=x_c(:,cell_neighbor(1));
%                 x_c_2=x_c(:,cell_neighbor(2));
%                 x_c_3=x_c(:,cell_neighbor(3));
%                 F=[(x_c_1-x_c(:,s_cell))';(x_c_2-x_c(:,s_cell))';(x_c_3-x_c(:,s_cell))'];
%             else
%                 error('logic error!');
%             end
% %             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';   
%         elseif FPDC==3 % All boundaries are periodic
%             cell_neighbor=cell_center{45};
%             CoeA=cell_center{48};
%             F=[(x_c(:,cell_neighbor(1))-x_c(:,s_cell))';(x_c(:,cell_neighbor(2))-x_c(:,s_cell))';(x_c(:,cell_neighbor(3))-x_c(:,s_cell))'];
% %             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             x=x_c(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeA*F))';
%         else
%             error('Logic error!');
%         end
%     else
%         error('Logic error!');
%     end
else
    error('The flag for mapping method is incorrect!');
end