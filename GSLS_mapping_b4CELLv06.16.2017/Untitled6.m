e=1000;
for o=1:O
    FC=FACE{o};
    if FC{10}==0 && FC{11}==0
        UP=FC{12};
        DOWN=FC{13};
        UP_Cell_number=UP(1,1);
        DOWN_Cell_number=DOWN(1,1);
        Cell_up_fixed=CELL{UP_Cell_number};
        C_up_fixed=Cell_up_fixed{5};
        Cell_down_fixed=CELL{DOWN_Cell_number};
        C_down_fixed=Cell_down_fixed{5};
        ND1=NODE{FC{8}};
        ND2=NODE{FC{9}};
        Cell_union=union(ND1{5},ND2{5});
        Cell_union_UP=setxor(UP_Cell_number,Cell_union);
        Cell_union_DOWN=setxor(DOWN_Cell_number,Cell_union);
        L_cell_union_up=length(Cell_union_UP);
        L_cell_union_down=length(Cell_union_DOWN);
        Pair_cell_union_up=zeros(2,L_cell_union_up*(L_cell_union_up-1)/2);
        Pair_cell_union_down=zeros(2,L_cell_union_down*(L_cell_union_down-1)/2);
        S=FC{18};
        C_down=S(:,1);
        C_up=S(:,2);
        %% Find the triangle that circles the stencil point in the upstream cell
        if single(e+norm(C_up_fixed-C_up))==single(e) % The stencil point is at the cnetroid
            ;
        else
            %% Find all possible paired centroids that could be used to form triangles
            a_up=0;
            for k=1:length(Cell_union_UP)
                if length(Cell_union_UP)==2
                    a_up=a_up+1;
                    Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(2)];
                    break;
                else
                    for l=2:length(Cell_union_UP)
                        a_up=a_up+1;
                        Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(l)];
                    end
                    Cell_union_UP=setxor(Cell_union_UP(1),Cell_union_UP);
                end
            end
            %% Narrow down the previous possibility ruling out the triangles that don't
            % circle the stencil point in the upstream cell
            a=0;
            in_tri_pair_up=0;
            for k=1:a_up
                Cell_pair_1=CELL{Pair_cell_union_up(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_up(2,k)};
                C_up_1=Cell_pair_1{5};
                C_up_2=Cell_pair_2{5};
                if in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                    a=a+1;
                    in_tri_pair_up(a)=k;
                end
            end
            if a==0
                error('No triangle is found that contains the stencil point!');
            end
            
            Pair_cell_union_up_new=zeros(2,a);
            for l=1:a
                Pair_cell_union_up_new(:,l)=Pair_cell_union_up(:,in_tri_pair_up(l));
            end
            % Check
            for k=1:a
                Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
                C_up_1=Cell_pair_1{5};
                C_up_2=Cell_pair_2{5};
                if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                    error(':');
                end
            end
            %% Find the pair that has the shortest distance to the stencil point
            Dis_up=zeros(1,a);
            for k=1:a
                Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
                C_up_1=Cell_pair_1{5};
                C_up_2=Cell_pair_2{5};
                Dis_up(1,k)=(dis(C_up,C_up_1)+dis(C_up,C_up_2))/2;
            end
            Dis_up_min=min(Dis_up);
            for k=1:a
                if single(e+Dis_up_min)==single(e+Dis_up(1,k))
                    break;
                end
            end
            Cell_pair_up_found=Pair_cell_union_up_new(:,k);
            % Check
            Cell_pair_1=CELL{Cell_pair_up_found(1,1)};
            Cell_pair_2=CELL{Cell_pair_up_found(2,1)};
            C_up_1=Cell_pair_1{5};
            C_up_2=Cell_pair_2{5};
            if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                error(':');
            end
        end
        
        %% Find the triangle that circles the stencil point in the downstream cell
        if single(e+norm(C_down_fixed-C_down))==single(e) % The stencil point is at the cnetroid
            ;
        else
            %% Find all possible paired centroids that could be used to form triangles
            a_down=0;
            for k=1:length(Cell_union_DOWN)
                if length(Cell_union_DOWN)==2
                    a_down=a_down+1;
                    Pair_cell_union_down(:,a_down)=[Cell_union_DOWN(1);Cell_union_DOWN(2)];
                    break;
                else
                    for l=2:length(Cell_union_DOWN)
                        a_down=a_down+1;
                        Pair_cell_union_down(:,a_down)=[Cell_union_DOWN(1);Cell_union_DOWN(l)];
                    end
                    Cell_union_DOWN=setxor(Cell_union_DOWN(1),Cell_union_DOWN);
                end
            end
            %% Narrow down the previous possibility ruling out the triangles that don't
            % circle the stencil point in the downstream cell
            a=0;
            in_tri_pair_down=0;
            for k=1:a_up
                Cell_pair_1=CELL{Pair_cell_union_down(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_down(2,k)};
                C_down_1=Cell_pair_1{5};
                C_down_2=Cell_pair_2{5};
                if in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                    a=a+1;
                    in_tri_pair_down(a)=k;
                end
            end
            if a==0
                error('No triangle is found that contains the stencil point!');
            end
            
            Pair_cell_union_down_new=zeros(2,a);
            for l=1:a
                Pair_cell_union_down_new(:,l)=Pair_cell_union_down(:,in_tri_pair_down(l));
            end
            % Check
            for k=1:a
                Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
                C_down_1=Cell_pair_1{5};
                C_down_2=Cell_pair_2{5};
                if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                    error(':');
                end
            end
            %% Find the pair that has the shortest distance to the stencil point
            Dis_down=zeros(1,a);
            for k=1:a
                Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
                Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
                C_down_1=Cell_pair_1{5};
                C_down_2=Cell_pair_2{5};
                Dis_down(1,k)=(dis(C_down,C_down_1)+dis(C_down,C_down_2))/2;
            end
            Dis_down_min=min(Dis_down);
            for k=1:a
                if single(e+Dis_down_min)==single(e+Dis_down(1,k))
                    break;
                end
            end
            Cell_pair_down_found=Pair_cell_union_down_new(:,k);
            % Check
            Cell_pair_1=CELL{Cell_pair_down_found(1,1)};
            Cell_pair_2=CELL{Cell_pair_down_found(2,1)};
            C_down_1=Cell_pair_1{5};
            C_down_2=Cell_pair_2{5};
            if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                error(':');
            end
        end
    end
end