function [result, current_cell_num_new] = point_value(Coord,current_cell_num_old,value_cell,value_nd,CELL,M,NODE,FPDC)
% result = point_value(Coord,X1,X2,Y1,Y2,value_cell,value_nd,CELL,M)
% returns the evaluated value at any given
% point in the computation domain.
% Coord is the coordinates of the given point, has to be column vector
% X1, X2, Y1 and Y2 is the boundaries of outer boundaries
% value_cell is the value from the centroid of triangles that is going to be
% sampled.
% value_nd is the value from the nodes that is going to be sampled.
% CELL is the cell data structure of all triangles
% M is the total number of triangles




P=CELL{current_cell_num_old};
if in_triangle(Coord,P{13},P{14},P{15})
%     [result,Cxy]=pfls(Coord,current_cell_num_old,CELL,value_cell,value_nd,FPDC);
    result=value_cell(:,current_cell_num_old);
    current_cell_num_new=current_cell_num_old;
else
    nd1=P{7};
    nd2=P{8};
    nd3=P{9};
    nd1_coord=P{13};
    nd2_coord=P{14};
    nd3_coord=P{15};
    if on_edge(Coord,nd1_coord,nd2_coord)
%         w1=dis(Coord,nd1_coord)/dis(nd1_coord,nd2_coord);
%         w2=1-w1;
%         result=w2*value_nd(:,nd1)+w1*value_nd(:,nd2);
        [result,Cxy]=pfls(Coord,current_cell_num_old,CELL,value_cell,value_nd,FPDC);

        ND1=NODE{nd1};
        ND2=NODE{nd2};
        cell_common=intersect(ND1{5},ND2{5});
        if cell_common==0
            error('Logic Error!');
        else
            current_cell_num_new=setxor(cell_common,current_cell_num_old);
            if length(current_cell_num_new)~=1
                error('Logic Error!');
            end
        end
    elseif on_edge(Coord,nd2_coord,nd3_coord)
%         w1=dis(Coord,nd2_coord)/dis(nd2_coord,nd3_coord);
%         w2=1-w1;
%         result=w2*value_nd(:,nd2)+w1*value_nd(:,nd3);
        [result,Cxy]=pfls(Coord,current_cell_num_old,CELL,value_cell,value_nd,FPDC);

        ND2=NODE{nd2};
        ND3=NODE{nd3};
        cell_common=intersect(ND2{5},ND3{5});
        if cell_common==0
            error('Logic Error!');
        else
            current_cell_num_new=setxor(cell_common,current_cell_num_old);
            if length(current_cell_num_new)~=1
                error('Logic Error!');
            end
        end
    elseif on_edge(Coord,nd3_coord,nd1_coord)
%         w1=dis(Coord,nd3_coord)/dis(nd3_coord,nd1_coord);
%         w2=1-w1;
%         result=w2*value_nd(:,nd3)+w1*value_nd(:,nd1);
        [result,Cxy]=pfls(Coord,current_cell_num_old,CELL,value_cell,value_nd,FPDC);

        ND3=NODE{nd3};
        ND1=NODE{nd1};
        cell_common=intersect(ND3{5},ND1{5});
        if cell_common==0
            error('Logic Error!');
        else
            current_cell_num_new=setxor(cell_common,current_cell_num_old);
            if length(current_cell_num_new)~=1
                error('Logic Error!');
            end
        end
    else
        cell_neigh=[P{2},P{3},P{4}];
        if length(setxor(cell_neigh,0))==4
            ;
        else
            cell_neigh=setxor(cell_neigh,0);
        end
        counter=0;
        for i=1:length(cell_neigh)
            Q=CELL{cell_neigh(i)};
            if in_triangle(Coord,Q{13},Q{14},Q{15})
                counter=counter+1;
                [result,Cxy]=pfls(Coord,cell_neigh(i),CELL,value_cell,value_nd,FPDC);
                current_cell_num_new=cell_neigh(i);
                break;
            end
        end
        if counter==1
            ;
        elseif counter==0
            r=in_which_triangle(Coord,CELL,M);
            if r==0 % The point is not within any triangle
                [nd1,nd2,ratio]=on_which_edge(Coord,CELL,M);
                if nd1==0 || nd2==0 % The given point is not on any edge
                    DIS=zeros(1,M);
                    for i=1:M
                        CL=CELL{i};
                        DIS(i)=dis(CL{5},Coord);
                    end
                    DIS_min=min(DIS);
                    for i=1:M
                        if single(DIS(i))==single(DIS_min)
                            break;
                        end
                    end
                    if i==M
                        if DIS(i)~=DIS_min
                            error('Logic Error!');
                        end
                    end
                    current_cell_num_new=i;
                    result=value_cell(:,i);
%                     [result,Cxy]=pfls(Coord,current_cell_num_new,CELL,value_cell,value_nd,FPDC);
%                     error('Logic Error!');
                else
%                     result=(1-ratio)*value_nd(:,nd1)+ratio*value_nd(:,nd2);
                    [result,Cxy]=pfls(Coord,current_cell_num_old,CELL,value_cell,value_nd,FPDC);
                    ND1=NODE{nd1};
                    ND2=NODE{nd2};
                    cell_common=intersect(ND1{5},ND2{5});
                    if cell_common==0
                        error('Logic Error!');
                    else
                        current_cell_num_new=cell_common(randi([1,2]));
                    end
                end
            else
                [result,Cxy]=pfls(Coord,r,CELL,value_cell,value_nd,FPDC);
                current_cell_num_new=r;
            end
        else
            error('Logic Error!');
        end
    end
end


% %% Old scheme, very slow
% e=((Y2-Y1)+(X2-X1))/2;
% e=2;
% e1=1e-7;
% r=in_which_triangle(Coord,CELL,M);
% if r==0 % The point is not within any triangle
%     [nd1,nd2,ratio]=on_which_edge(Coord,CELL,M);
%     if nd1==0 || nd2==0 % The given point is not on any edge
%         disp('The given point is out of computational domain!');
%         result=0;
%     else
%         result=(1-ratio)*value_nd(:,nd1)+ratio*value_nd(:,nd2);
%     end
% else % Found the triangle that enclose the point
%     P=CELL{r};
%     if dis(P{5},Coord)<e1
%         in_zone_four=1;
%     else
%         in_zone_one=in_triangle(Coord,P{5},P{13},P{14});
%         in_zone_two=in_triangle(Coord,P{5},P{14},P{15});
%         in_zone_three=in_triangle(Coord,P{5},P{15},P{13});
%         in_zone_four=0;
%     end
%     
%     if in_zone_four
%         result=value_cell(:,P{1}); % follow the assumption that the value in the triangle is constant
%     else
%         if in_zone_one
%             result=in_cell_mapping_zone(P,1,value_cell(:,P{1}),value_nd(:,P{7}),value_nd(:,P{8}),value_nd(:,P{9}),Coord,2); % Located in zone 1
%         elseif in_zone_two
%             result=in_cell_mapping_zone(P,2,value_cell(:,P{1}),value_nd(:,P{7}),value_nd(:,P{8}),value_nd(:,P{9}),Coord,2); % Located in zone 2
%         elseif in_zone_three
%             result=in_cell_mapping_zone(P,3,value_cell(:,P{1}),value_nd(:,P{7}),value_nd(:,P{8}),value_nd(:,P{9}),Coord,2); % Located in zone 3
%         else
%             on_edge_one=on_edge(Coord,P{5},P{13});
%             on_edge_two=on_edge(Coord,P{5},P{14});
%             on_edge_three=on_edge(Coord,P{5},P{15});
%             % Check
%             if (on_edge_one+on_edge_two+on_edge_three)~=1 && (on_edge_one+on_edge_two+on_edge_three)~=3
%                 error('The stencil point could no be on multiple edges!');
%             end
%             if on_edge_one && (on_edge_two && on_edge_three)
%                 % Check
%                 if single(e+dis(Coord,P{5}))~=single(e)
%                     error('Logic error!');
%                 end
%                 result=value_cell(:,P{1}); % follow the assumption that the value in the triangle is constant
%             elseif on_edge_one
%                 ratio=dis(Coord,P{5})/dis(P{5},P{13});
%                 result=(1-ratio)*value_cell(:,P{1})+ratio*value_nd(:,P{7});
%             elseif on_edge_two
%                 ratio=dis(Coord,P{5})/dis(P{5},P{14});
%                 result=(1-ratio)*value_cell(:,P{1})+ratio*value_nd(:,P{8});
%             elseif on_edge_three
%                 ratio=dis(Coord,P{5})/dis(P{5},P{15});
%                 result=(1-ratio)*value_cell(:,P{1})+ratio*value_nd(:,P{9});
%             else
%                 error('The zone for the stencil point is not found!');
%             end
%         end
%     end
% end