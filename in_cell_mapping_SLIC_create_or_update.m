function [CELL,NODE,FACE]=in_cell_mapping_SLIC_create_or_update(CELL,NODE,FACE,X1,X2,Y1,Y2,V,dt)

%% This function creat the pre-computed information for the application of second-order interpolation scheme for the Semi-Lagrangian Implicit Collision (SLIC) scheme.

M=length(CELL);
[d,q]=size(V);
SL_location_coord=zeros(d,q);

%% Algorithm 1&2, equally weighted 2nd-order interpolation see CELLv.11.11.2019
%% Fill CL{57~58} for cell type 1 and CL{59~61} for cell type 2 & 3
for i=1:M
    CL=CELL{i};
    if CL{37}==1 % cell type 1
        cell1=CL{2};
        cell2=CL{3};
        cell3=CL{4};
        Cell1=CELL{cell1};
        Cell2=CELL{cell2};
        Cell3=CELL{cell3};
        
        A=[Cell1{5}-CL{5},Cell2{5}-CL{5},Cell3{5}-CL{5}]';
        
        for k=1:q
            SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
            if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
                error('Please decrease the size of dt!');
            end
        end
        
        CL{57}=[CL{2},CL{3},CL{4}];
        CL{58}=((SL_location_coord-CL{5})')*(inv(A'*A)*A');
    else % cell type 2 and 3
        cell1=CL{2};
        cell2=CL{3};
        cell3=CL{4};
        %% Usung three point, one of which is from boundary
        if cell1==0
            FC_bc=FACE{CL{16}};
            Cell2=CELL{cell2};
            Cell3=CELL{cell3};
                            coor_bc=CL{28};
                            coor_1=coor_bc;
%             coor_1=mirror(CL{5},FC_bc{5},FC_bc{6});
            coor_2=Cell2{5};
            coor_3=Cell3{5};
        elseif cell2==0
            FC_bc=FACE{CL{17}};
            Cell1=CELL{cell1};
            Cell3=CELL{cell3};
                            coor_bc=CL{29};
            coor_1=Cell1{5};
                            coor_2=coor_bc;
%             coor_2=mirror(CL{5},FC_bc{5},FC_bc{6});
            coor_3=Cell3{5};
        elseif cell3==0
            FC_bc=FACE{CL{18}};
            Cell1=CELL{cell1};
            Cell2=CELL{cell2};
                            coor_bc=CL{30};
            coor_1=Cell1{5};
            coor_2=Cell2{5};
                            coor_3=coor_bc;
%             coor_3=mirror(CL{5},FC_bc{5},FC_bc{6});
        else
            error('Logic error!');
        end
        
        %             A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]'/sqrt(CL{6});
        A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]';
        %% Using only existing centroids, don't use boundary node at all
        %             if cell1==0
        %                 FC_bc=FACE{CL{16}};
        %                 Cell2=CELL{cell2};
        %                 Cell3=CELL{cell3};
        %                 coor_2=Cell2{5};
        %                 coor_3=Cell3{5};
        %                 A=[coor_2-CL{5},coor_3-CL{5}]';
        %             elseif cell2==0
        %                 FC_bc=FACE{CL{17}};
        %                 Cell1=CELL{cell1};
        %                 Cell3=CELL{cell3};
        %                 coor_1=Cell1{5};
        %                 coor_3=Cell3{5};
        %                 A=[coor_1-CL{5},coor_3-CL{5}]';
        %             elseif cell3==0
        %                 FC_bc=FACE{CL{18}};
        %                 Cell1=CELL{cell1};
        %                 Cell2=CELL{cell2};
        %                 coor_1=Cell1{5};
        %                 coor_2=Cell2{5};
        %                 A=[coor_1-CL{5},coor_2-CL{5}]';
        %             else
        %                 error('Logic error!');
        %             end
        
        
        CL{59}=[CL{2},CL{3},CL{4}];
        
        bc_intercept=FC_bc{18}; % See definition in FACE v03.12.2016
        bc_intcpt_1=bc_intercept(:,2);
        bc_intcpt_2=bc_intercept(:,3);
        if sum(bc_intcpt_1)==0 && sum(bc_intcpt_2)~=0
            if bc_intcpt_2(1)==FC_bc{8} && bc_intcpt_2(2)==FC_bc{9}
                ;
            elseif bc_intcpt_2(1)==FC_bc{9} && bc_intcpt_2(2)==FC_bc{8}
                ;
            else
                error('Logic error!');
            end
            CL{60}=bc_intcpt_2';
        elseif sum(bc_intcpt_1)~=0 && sum(bc_intcpt_2)==0
            if bc_intcpt_1(1)==FC_bc{8} && bc_intcpt_1(2)==FC_bc{9}
                ;
            elseif bc_intcpt_1(1)==FC_bc{9} && bc_intcpt_1(2)==FC_bc{8}
                ;
            else
                error('Logic error!');
            end
            CL{60}=bc_intcpt_1';
        else
            error('Logic Error!');
        end
        
        for k=1:q
            SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
            if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
                error('Please decrease the size of dt!');
            end
        end
        
        CL{61}=((SL_location_coord-CL{5})')*(inv(A'*A)*A');
    end
    CELL{i}=CL;
end


%% Fill CL{62~64}
for i=1:M
    CL=CELL{i};
    if CL{37}==2
        if CL{2}==0
            FC_on_bc=FACE{CL{16}};
        elseif CL{3}==0
            FC_on_bc=FACE{CL{17}};
        elseif CL{4}==0
            FC_on_bc=FACE{CL{18}};
        else
            error('Logic error');
        end
        if FC_on_bc{24}==2 % The face on outer boundary
            centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
            if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
                error('logic error');
            else
                cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
            end
            bc_cross_stencil_line=FC_on_bc{23};
            if length(setxor(bc_cross_stencil_line(2:3),0))~=1
                error('logic error');
            else
                bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
            end
            if bc_cross_id==1
                centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
            elseif bc_cross_id==2
                centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
            elseif bc_cross_id==3
                centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
            elseif bc_cross_id==4
                centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
            else
                error('logic error!');
            end
            
            if CL{2}==0
                neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
                cell2=CL{3};
                cell3=CL{4};
                Cell2=CELL{cell2};
                Cell3=CELL{cell3};
                Coor1=centroid_cell_neigh_on_bc;
                Coor2=Cell2{5};
                Coor3=Cell3{5};
            elseif CL{3}==0
                neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
                cell1=CL{2};
                cell3=CL{4};
                Cell1=CELL{cell1};
                Cell3=CELL{cell3};
                Coor1=Cell1{5};
                Coor2=centroid_cell_neigh_on_bc;
                Coor3=Cell3{5};
            elseif CL{4}==0
                neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
                cell1=CL{2};
                cell2=CL{3};
                Cell1=CELL{cell1};
                Cell2=CELL{cell2};
                Coor1=Cell1{5};
                Coor2=Cell2{5};
                Coor3=centroid_cell_neigh_on_bc;
            else
                error('Logic error');
            end
            
            %                 A=[Coor1-CL{5},Coor2-CL{5},Coor3-CL{5}]'/sqrt(CL{6});
            A=[Coor1-CL{5},Coor2-CL{5},Coor3-CL{5}]';
            
            for k=1:q
                SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
                if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
                    error('Please decrease the size of dt!');
                end
            end
            
            CL{62}=neigh_id;
            CL{63}=[0,0,0];
            CL{64}=((SL_location_coord-CL{5})')*(inv(A'*A)*A');
            
            CELL{i}=CL;
        else
            error('logic error!');
        end
    end
end



% %% Algorithm 3&4, upwind 2nd-order interpolation
% %% Fill CL{57~59} for cell type 1 and CL{60~63} for cell type 2 & 3
% for i=1:M
%     CL=CELL{i};
%     if CL{37}==1 % cell type 1
%         cell1=CL{2};
%         cell2=CL{3};
%         cell3=CL{4};
%         Cell1=CELL{cell1};
%         Cell2=CELL{cell2};
%         Cell3=CELL{cell3};
%         
%         A1=[(Cell1{5}-CL{5}),Cell2{5}-CL{5}]';
%         A2=[(Cell2{5}-CL{5}),Cell3{5}-CL{5}]';
%         A3=[(Cell3{5}-CL{5}),Cell1{5}-CL{5}]';
%         
%         Zone_ID=zeros(1,q);
%         
%         A_matrix=zeros(q,2);
%         DIS_01=dis(CL{5},Cell1{5});
%         DIS_02=dis(CL{5},Cell2{5});
%         DIS_03=dis(CL{5},Cell3{5});
%         e=(DIS_01+DIS_02+DIS_03)/3*10;
%         for k=1:q
%             in_zone_one=0;
%             in_zone_two=0;
%             in_zone_three=0;
%             SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
%             if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
%                 error('Please decrease the size of dt!');
%             else
%                 in_zone_one=in_triangle(SL_location_coord(:,k),CL{5},Cell1{5},Cell2{5});
%                 in_zone_two=in_triangle(SL_location_coord(:,k),CL{5},Cell2{5},Cell3{5});
%                 in_zone_three=in_triangle(SL_location_coord(:,k),CL{5},Cell3{5},Cell1{5});
%                 if (in_zone_one+in_zone_two+in_zone_three)==1
%                     if in_zone_one==1
%                         Zone_ID(1,k)=1;
%                         A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A1'*A1)*A1');
%                     elseif in_zone_two==1
%                         Zone_ID(1,k)=2;
%                         A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A2'*A2)*A2');
%                     elseif in_zone_three==1
%                         Zone_ID(1,k)=3;
%                         A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A3'*A3)*A3');
%                     else
%                         error('Logic Error!');
%                     end
%                 else
%                     if single(dis(SL_location_coord(:,k),CL{5})+e)==single(e)
%                         if k==1
%                             Zone_ID(1,k)=123;
%                             A_matrix(k,:)=[0,0];
%                         else
%                             error('Logic Error!');
%                         end
%                     else
%                         if single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell1{5})+e)==single(DIS_01+e)
%                             Zone_ID(1,k)=31;
%                             A_matrix(k,:)=[dis(SL_location_coord(:,k),Cell1{5}),dis(SL_location_coord(:,k),CL{5})]/DIS_01;
%                         elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell2{5})+e)==single(DIS_02+e)
%                             Zone_ID(1,k)=12;
%                             A_matrix(k,:)=[dis(SL_location_coord(:,k),Cell2{5}),dis(SL_location_coord(:,k),CL{5})]/DIS_02;
%                         elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell3{5})+e)==single(DIS_03+e)
%                             Zone_ID(1,k)=23;
%                             A_matrix(k,:)=[dis(SL_location_coord(:,k),Cell3{5}),dis(SL_location_coord(:,k),CL{5})]/DIS_03;
%                         else
%                             error('Logic Error!');
%                         end
%                     end
%                 end
%             end
%         end
%         
%         CL{57}=[CL{2},CL{3},CL{4}];
%         CL{58}=A_matrix;
%         CL{59}=Zone_ID;
%     else % cell type 2 and 3
% %         cell1=CL{2};
% %         cell2=CL{3};
% %         cell3=CL{4};
% %         %% Usung three point, one of which is from boundary
% %         if cell1==0
% %             FC_bc=FACE{CL{16}};
% %             Cell2=CELL{cell2};
% %             Cell3=CELL{cell3};
% %             %                 coor_bc=CL{28};
% %             %                 coor_1=coor_bc;
% %             coor_1=mirror(CL{5},FC_bc{5},FC_bc{6});
% %             coor_2=Cell2{5};
% %             coor_3=Cell3{5};
% %         elseif cell2==0
% %             FC_bc=FACE{CL{17}};
% %             Cell1=CELL{cell1};
% %             Cell3=CELL{cell3};
% %             %                 coor_bc=CL{29};
% %             coor_1=Cell1{5};
% %             %                 coor_2=coor_bc;
% %             coor_2=mirror(CL{5},FC_bc{5},FC_bc{6});
% %             coor_3=Cell3{5};
% %         elseif cell3==0
% %             FC_bc=FACE{CL{18}};
% %             Cell1=CELL{cell1};
% %             Cell2=CELL{cell2};
% %             %                 coor_bc=CL{30};
% %             coor_1=Cell1{5};
% %             coor_2=Cell2{5};
% %             %                 coor_3=coor_bc;
% %             coor_3=mirror(CL{5},FC_bc{5},FC_bc{6});
% %         else
% %             error('Logic error!');
% %         end
% %         
% %         %             A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]'/sqrt(CL{6});
% %         A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]';
% %         %% Using only existing centroids, don't use boundary node at all
% %         %             if cell1==0
% %         %                 FC_bc=FACE{CL{16}};
% %         %                 Cell2=CELL{cell2};
% %         %                 Cell3=CELL{cell3};
% %         %                 coor_2=Cell2{5};
% %         %                 coor_3=Cell3{5};
% %         %                 A=[coor_2-CL{5},coor_3-CL{5}]';
% %         %             elseif cell2==0
% %         %                 FC_bc=FACE{CL{17}};
% %         %                 Cell1=CELL{cell1};
% %         %                 Cell3=CELL{cell3};
% %         %                 coor_1=Cell1{5};
% %         %                 coor_3=Cell3{5};
% %         %                 A=[coor_1-CL{5},coor_3-CL{5}]';
% %         %             elseif cell3==0
% %         %                 FC_bc=FACE{CL{18}};
% %         %                 Cell1=CELL{cell1};
% %         %                 Cell2=CELL{cell2};
% %         %                 coor_1=Cell1{5};
% %         %                 coor_2=Cell2{5};
% %         %                 A=[coor_1-CL{5},coor_2-CL{5}]';
% %         %             else
% %         %                 error('Logic error!');
% %         %             end
% %         
% %         
% %         CL{59}=[CL{2},CL{3},CL{4}];
% %         
% %         bc_intercept=FC_bc{18}; % See definition in FACE v03.12.2016
% %         bc_intcpt_1=bc_intercept(:,2);
% %         bc_intcpt_2=bc_intercept(:,3);
% %         if sum(bc_intcpt_1)==0 && sum(bc_intcpt_2)~=0
% %             if bc_intcpt_2(1)==FC_bc{8} && bc_intcpt_2(2)==FC_bc{9}
% %                 ;
% %             elseif bc_intcpt_2(1)==FC_bc{9} && bc_intcpt_2(2)==FC_bc{8}
% %                 ;
% %             else
% %                 error('Logic error!');
% %             end
% %             CL{60}=bc_intcpt_2';
% %         elseif sum(bc_intcpt_1)~=0 && sum(bc_intcpt_2)==0
% %             if bc_intcpt_1(1)==FC_bc{8} && bc_intcpt_1(2)==FC_bc{9}
% %                 ;
% %             elseif bc_intcpt_1(1)==FC_bc{9} && bc_intcpt_1(2)==FC_bc{8}
% %                 ;
% %             else
% %                 error('Logic error!');
% %             end
% %             CL{60}=bc_intcpt_1';
% %         else
% %             error('Logic Error!');
% %         end
% %         
% %         for k=1:q
% %             SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
% %             if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
% %                 error('Please decrease the size of dt!');
% %             end
% %         end
% %         
% %         CL{61}=((SL_location_coord-CL{5})')*(inv(A'*A)*A');
%     end
%     CELL{i}=CL;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %% Fill CL{64~67}
% for i=1:M
%     CL=CELL{i};
%     if CL{37}==2
%         if CL{2}==0
%             FC_on_bc=FACE{CL{16}};
%         elseif CL{3}==0
%             FC_on_bc=FACE{CL{17}};
%         elseif CL{4}==0
%             FC_on_bc=FACE{CL{18}};
%         else
%             error('Logic error');
%         end
%         if FC_on_bc{24}==2 % The face on outer boundary
%             centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
%             if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
%                 error('logic error');
%             else
%                 cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
%             end
%             bc_cross_stencil_line=FC_on_bc{23};
%             if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                 error('logic error');
%             else
%                 bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%             end
%             if bc_cross_id==1
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
%             elseif bc_cross_id==2
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
%             elseif bc_cross_id==3
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
%             elseif bc_cross_id==4
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
%             else
%                 error('logic error!');
%             end
%             
%             if CL{2}==0
%                 neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
%                 cell2=CL{3};
%                 cell3=CL{4};
%                 Cell2=CELL{cell2};
%                 Cell3=CELL{cell3};
%                 Coor1=centroid_cell_neigh_on_bc;
%                 Coor2=Cell2{5};
%                 Coor3=Cell3{5};
%             elseif CL{3}==0
%                 neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
%                 cell1=CL{2};
%                 cell3=CL{4};
%                 Cell1=CELL{cell1};
%                 Cell3=CELL{cell3};
%                 Coor1=Cell1{5};
%                 Coor2=centroid_cell_neigh_on_bc;
%                 Coor3=Cell3{5};
%             elseif CL{4}==0
%                 neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
%                 cell1=CL{2};
%                 cell2=CL{3};
%                 Cell1=CELL{cell1};
%                 Cell2=CELL{cell2};
%                 Coor1=Cell1{5};
%                 Coor2=Cell2{5};
%                 Coor3=centroid_cell_neigh_on_bc;
%             else
%                 error('Logic error');
%             end
%             
%             %                 A=[Coor1-CL{5},Coor2-CL{5},Coor3-CL{5}]'/sqrt(CL{6});
%             
%             A1=[(Coor1-CL{5}),Coor2-CL{5}]';
%             A2=[(Coor2-CL{5}),Coor3-CL{5}]';
%             A3=[(Coor3-CL{5}),Coor1-CL{5}]';
%             
%             Zone_ID=zeros(1,q);
%             
%             A_matrix=zeros(q,2);
%             DIS_01=dis(CL{5},Coor1);
%             DIS_02=dis(CL{5},Coor2);
%             DIS_03=dis(CL{5},Coor3);
%             e=(DIS_01+DIS_02+DIS_03)/3*10;
%             for k=1:q
%                 in_zone_one=0;
%                 in_zone_two=0;
%                 in_zone_three=0;
%                 SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
%                 if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
%                     error('Please decrease the size of dt!');
%                 else
%                     in_zone_one=in_triangle(SL_location_coord(:,k),CL{5},Coor1,Coor2);
%                     in_zone_two=in_triangle(SL_location_coord(:,k),CL{5},Coor2,Coor3);
%                     in_zone_three=in_triangle(SL_location_coord(:,k),CL{5},Coor3,Coor1);
%                     if (in_zone_one+in_zone_two+in_zone_three)==1
%                         if in_zone_one==1
%                             Zone_ID(1,k)=1;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A1'*A1)*A1');
%                         elseif in_zone_two==1
%                             Zone_ID(1,k)=2;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A2'*A2)*A2');
%                         elseif in_zone_three==1
%                             Zone_ID(1,k)=3;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A3'*A3)*A3');
%                         else
%                             error('Logic Error!');
%                         end
%                     else
%                         if single(dis(SL_location_coord(:,k),CL{5})+e)==single(e)
%                             if k==1
%                                 Zone_ID(1,k)=123;
%                                 A_matrix(k,:)=[0,0];
%                             else
%                                 error('Logic Error!');
%                             end
%                         else
%                             if single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor1)+e)==single(DIS_01+e)
%                                 Zone_ID(1,k)=31;
%                                 A_matrix(k,:)=[dis(SL_location_coord(:,k),Coor1),dis(SL_location_coord(:,k),CL{5})]/DIS_01;
%                             elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor2)+e)==single(DIS_02+e)
%                                 Zone_ID(1,k)=12;
%                                 A_matrix(k,:)=[dis(SL_location_coord(:,k),Coor2),dis(SL_location_coord(:,k),CL{5})]/DIS_02;
%                             elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor3)+e)==single(DIS_03+e)
%                                 Zone_ID(1,k)=23;
%                                 A_matrix(k,:)=[dis(SL_location_coord(:,k),Coor3),dis(SL_location_coord(:,k),CL{5})]/DIS_03;
%                             else
%                                 error('Logic Error!');
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             CL{64}=neigh_id;
%             CL{65}=[0,0,0];
%             CL{66}=A_matrix;
%             CL{67}=Zone_ID;
%             
%             CELL{i}=CL;
%         else
%             error('logic error!');
%         end
%     end
% end


% %% Algorithm 5, downwind 2nd-order interpolation
% %% Fill CL{57~59} for cell type 1 and CL{60~63} for cell type 2 & 3
% for i=1:M
%     CL=CELL{i};
%     if CL{37}==1 % cell type 1
%         cell1=CL{2};
%         cell2=CL{3};
%         cell3=CL{4};
%         Cell1=CELL{cell1};
%         Cell2=CELL{cell2};
%         Cell3=CELL{cell3};
%         
%         A1=[(Cell1{5}-CL{5}),Cell2{5}-CL{5}]';
%         A2=[(Cell2{5}-CL{5}),Cell3{5}-CL{5}]';
%         A3=[(Cell3{5}-CL{5}),Cell1{5}-CL{5}]';
%         
%         Zone_ID=zeros(1,q);
%         
%         A_matrix=zeros(q,2);
%         DIS_01=dis(CL{5},Cell1{5});
%         DIS_02=dis(CL{5},Cell2{5});
%         DIS_03=dis(CL{5},Cell3{5});
%         e=(DIS_01+DIS_02+DIS_03)/3*10;
%         for k=1:q
%             in_zone_one=0;
%             in_zone_two=0;
%             in_zone_three=0;
%             SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
%             if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
%                 error('Please decrease the size of dt!');
%             else
%                 in_zone_one=in_triangle(SL_location_coord(:,k),CL{5},Cell1{5},Cell2{5});
%                 in_zone_two=in_triangle(SL_location_coord(:,k),CL{5},Cell2{5},Cell3{5});
%                 in_zone_three=in_triangle(SL_location_coord(:,k),CL{5},Cell3{5},Cell1{5});
%                 if (in_zone_one+in_zone_two+in_zone_three)==1
%                     if in_zone_one==1
%                         Zone_ID(1,k)=1;
%                         C_j=norm_joint(CL{5},CL{5}-Cell3{5},SL_location_coord(:,k));
%                         DIS_ex=dis(C_j,CL{5});
%                         if single(e+dis(C_j,Cell3{5}))~=single(e+DIS_ex+DIS_03)
%                             error('Logic Error!');
%                         else
%                             A_matrix(k,:)=[(DIS_ex+DIS_03)/DIS_03,-DIS_ex/DIS_03];
%                         end
%                     elseif in_zone_two==1
%                         Zone_ID(1,k)=2;
%                         C_j=norm_joint(CL{5},CL{5}-Cell1{5},SL_location_coord(:,k));
%                         DIS_ex=dis(C_j,CL{5});
%                         if single(e+dis(C_j,Cell1{5}))~=single(e+DIS_ex+DIS_01)
%                             error('Logic Error!');
%                         else
%                             A_matrix(k,:)=[(DIS_ex+DIS_01)/DIS_01,-DIS_ex/DIS_01];
%                         end
%                     elseif in_zone_three==1
%                         Zone_ID(1,k)=3;
%                         C_j=norm_joint(CL{5},CL{5}-Cell2{5},SL_location_coord(:,k));
%                         DIS_ex=dis(C_j,CL{5});
%                         if single(e+dis(C_j,Cell2{5}))~=single(e+DIS_ex+DIS_02)
%                             error('Logic Error!');
%                         else
%                             A_matrix(k,:)=[(DIS_ex+DIS_02)/DIS_02,-DIS_ex/DIS_02];
%                         end
%                     else
%                         error('Logic Error!');
%                     end
%                 else
%                     if single(dis(SL_location_coord(:,k),CL{5})+e)==single(e)
%                         if k==1
%                             Zone_ID(1,k)=123;
%                             A_matrix(k,:)=[0,0];
%                         else
%                             error('Logic Error!');
%                         end
%                     else
%                         if single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell1{5})+e)==single(DIS_01+e)
%                             Zone_ID(1,k)=31;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A2'*A2)*A2');
%                         elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell2{5})+e)==single(DIS_02+e)
%                             Zone_ID(1,k)=12;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A3'*A3)*A3');
%                         elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Cell3{5})+e)==single(DIS_03+e)
%                             Zone_ID(1,k)=23;
%                             A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A1'*A1)*A1');
%                         else
%                             error('Logic Error!');
%                         end
%                     end
%                 end
%             end
%         end
%         
%         CL{57}=[CL{2},CL{3},CL{4}];
%         CL{58}=A_matrix;
%         CL{59}=Zone_ID;
%     else % cell type 2 and 3
% %         cell1=CL{2};
% %         cell2=CL{3};
% %         cell3=CL{4};
% %         %% Usung three point, one of which is from boundary
% %         if cell1==0
% %             FC_bc=FACE{CL{16}};
% %             Cell2=CELL{cell2};
% %             Cell3=CELL{cell3};
% %             %                 coor_bc=CL{28};
% %             %                 coor_1=coor_bc;
% %             coor_1=mirror(CL{5},FC_bc{5},FC_bc{6});
% %             coor_2=Cell2{5};
% %             coor_3=Cell3{5};
% %         elseif cell2==0
% %             FC_bc=FACE{CL{17}};
% %             Cell1=CELL{cell1};
% %             Cell3=CELL{cell3};
% %             %                 coor_bc=CL{29};
% %             coor_1=Cell1{5};
% %             %                 coor_2=coor_bc;
% %             coor_2=mirror(CL{5},FC_bc{5},FC_bc{6});
% %             coor_3=Cell3{5};
% %         elseif cell3==0
% %             FC_bc=FACE{CL{18}};
% %             Cell1=CELL{cell1};
% %             Cell2=CELL{cell2};
% %             %                 coor_bc=CL{30};
% %             coor_1=Cell1{5};
% %             coor_2=Cell2{5};
% %             %                 coor_3=coor_bc;
% %             coor_3=mirror(CL{5},FC_bc{5},FC_bc{6});
% %         else
% %             error('Logic error!');
% %         end
% %         
% %         %             A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]'/sqrt(CL{6});
% %         A=[coor_1-CL{5},coor_2-CL{5},coor_3-CL{5}]';
% %         %% Using only existing centroids, don't use boundary node at all
% %         %             if cell1==0
% %         %                 FC_bc=FACE{CL{16}};
% %         %                 Cell2=CELL{cell2};
% %         %                 Cell3=CELL{cell3};
% %         %                 coor_2=Cell2{5};
% %         %                 coor_3=Cell3{5};
% %         %                 A=[coor_2-CL{5},coor_3-CL{5}]';
% %         %             elseif cell2==0
% %         %                 FC_bc=FACE{CL{17}};
% %         %                 Cell1=CELL{cell1};
% %         %                 Cell3=CELL{cell3};
% %         %                 coor_1=Cell1{5};
% %         %                 coor_3=Cell3{5};
% %         %                 A=[coor_1-CL{5},coor_3-CL{5}]';
% %         %             elseif cell3==0
% %         %                 FC_bc=FACE{CL{18}};
% %         %                 Cell1=CELL{cell1};
% %         %                 Cell2=CELL{cell2};
% %         %                 coor_1=Cell1{5};
% %         %                 coor_2=Cell2{5};
% %         %                 A=[coor_1-CL{5},coor_2-CL{5}]';
% %         %             else
% %         %                 error('Logic error!');
% %         %             end
% %         
% %         
% %         CL{59}=[CL{2},CL{3},CL{4}];
% %         
% %         bc_intercept=FC_bc{18}; % See definition in FACE v03.12.2016
% %         bc_intcpt_1=bc_intercept(:,2);
% %         bc_intcpt_2=bc_intercept(:,3);
% %         if sum(bc_intcpt_1)==0 && sum(bc_intcpt_2)~=0
% %             if bc_intcpt_2(1)==FC_bc{8} && bc_intcpt_2(2)==FC_bc{9}
% %                 ;
% %             elseif bc_intcpt_2(1)==FC_bc{9} && bc_intcpt_2(2)==FC_bc{8}
% %                 ;
% %             else
% %                 error('Logic error!');
% %             end
% %             CL{60}=bc_intcpt_2';
% %         elseif sum(bc_intcpt_1)~=0 && sum(bc_intcpt_2)==0
% %             if bc_intcpt_1(1)==FC_bc{8} && bc_intcpt_1(2)==FC_bc{9}
% %                 ;
% %             elseif bc_intcpt_1(1)==FC_bc{9} && bc_intcpt_1(2)==FC_bc{8}
% %                 ;
% %             else
% %                 error('Logic error!');
% %             end
% %             CL{60}=bc_intcpt_1';
% %         else
% %             error('Logic Error!');
% %         end
% %         
% %         for k=1:q
% %             SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
% %             if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
% %                 error('Please decrease the size of dt!');
% %             end
% %         end
% %         
% %         CL{61}=((SL_location_coord-CL{5})')*(inv(A'*A)*A');
%     end
%     CELL{i}=CL;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %% Fill CL{64~67}
% for i=1:M
%     CL=CELL{i};
%     if CL{37}==2
%         if CL{2}==0
%             FC_on_bc=FACE{CL{16}};
%         elseif CL{3}==0
%             FC_on_bc=FACE{CL{17}};
%         elseif CL{4}==0
%             FC_on_bc=FACE{CL{18}};
%         else
%             error('Logic error');
%         end
%         if FC_on_bc{24}==2 % The face on outer boundary
%             centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
%             if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
%                 error('logic error');
%             else
%                 cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
%             end
%             bc_cross_stencil_line=FC_on_bc{23};
%             if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                 error('logic error');
%             else
%                 bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%             end
%             if bc_cross_id==1
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
%             elseif bc_cross_id==2
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
%             elseif bc_cross_id==3
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
%             elseif bc_cross_id==4
%                 centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
%             else
%                 error('logic error!');
%             end
%             
%             if CL{2}==0
%                 neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
%                 cell2=CL{3};
%                 cell3=CL{4};
%                 Cell2=CELL{cell2};
%                 Cell3=CELL{cell3};
%                 Coor1=centroid_cell_neigh_on_bc;
%                 Coor2=Cell2{5};
%                 Coor3=Cell3{5};
%             elseif CL{3}==0
%                 neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
%                 cell1=CL{2};
%                 cell3=CL{4};
%                 Cell1=CELL{cell1};
%                 Cell3=CELL{cell3};
%                 Coor1=Cell1{5};
%                 Coor2=centroid_cell_neigh_on_bc;
%                 Coor3=Cell3{5};
%             elseif CL{4}==0
%                 neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
%                 cell1=CL{2};
%                 cell2=CL{3};
%                 Cell1=CELL{cell1};
%                 Cell2=CELL{cell2};
%                 Coor1=Cell1{5};
%                 Coor2=Cell2{5};
%                 Coor3=centroid_cell_neigh_on_bc;
%             else
%                 error('Logic error');
%             end
%             
%             %                 A=[Coor1-CL{5},Coor2-CL{5},Coor3-CL{5}]'/sqrt(CL{6});
%             
%             A1=[(Coor1-CL{5}),Coor2-CL{5}]';
%             A2=[(Coor2-CL{5}),Coor3-CL{5}]';
%             A3=[(Coor3-CL{5}),Coor1-CL{5}]';
%             
%             Zone_ID=zeros(1,q);
%             
%             A_matrix=zeros(q,2);
%             DIS_01=dis(CL{5},Coor1);
%             DIS_02=dis(CL{5},Coor2);
%             DIS_03=dis(CL{5},Coor3);
%             e=(DIS_01+DIS_02+DIS_03)/3*10;
%             for k=1:q
%                 in_zone_one=0;
%                 in_zone_two=0;
%                 in_zone_three=0;
%                 SL_location_coord(:,k)=CL{5}+(-V(:,k)*(dt));
%                 if ~in_triangle(SL_location_coord(:,k),CL{13},CL{14},CL{15})
%                     error('Please decrease the size of dt!');
%                 else
%                     in_zone_one=in_triangle(SL_location_coord(:,k),CL{5},Coor1,Coor2);
%                     in_zone_two=in_triangle(SL_location_coord(:,k),CL{5},Coor2,Coor3);
%                     in_zone_three=in_triangle(SL_location_coord(:,k),CL{5},Coor3,Coor1);
%                     if (in_zone_one+in_zone_two+in_zone_three)==1
%                         if in_zone_one==1
%                             Zone_ID(1,k)=1;
%                             C_j=norm_joint(CL{5},norm_v(CL{5}-Coor3),SL_location_coord(:,k));
%                             DIS_ex=dis(C_j,CL{5});
%                             if single(e+dis(C_j,Coor3))~=single(e+DIS_ex+DIS_03)
%                                 error('Logic Error!');
%                             else
%                                 A_matrix(k,:)=[(DIS_ex+DIS_03)/DIS_03,-DIS_ex/DIS_03];
%                             end
%                         elseif in_zone_two==1
%                             Zone_ID(1,k)=2;
%                             C_j=norm_joint(CL{5},norm_v(CL{5}-Coor1),SL_location_coord(:,k));
%                             DIS_ex=dis(C_j,CL{5});
%                             if single(e+dis(C_j,Coor1))~=single(e+DIS_ex+DIS_01)
%                                 error('Logic Error!');
%                             else
%                                 A_matrix(k,:)=[(DIS_ex+DIS_01)/DIS_01,-DIS_ex/DIS_01];
%                             end
%                         elseif in_zone_three==1
%                             Zone_ID(1,k)=3;
%                             C_j=norm_joint(CL{5},norm_v(CL{5}-Coor2),SL_location_coord(:,k));
%                             DIS_ex=dis(C_j,CL{5});
%                             if single(e+dis(C_j,Coor2))~=single(e+DIS_ex+DIS_02)
%                                 error('Logic Error!');
%                             else
%                                 A_matrix(k,:)=[(DIS_ex+DIS_02)/DIS_02,-DIS_ex/DIS_02];
%                             end
%                         else
%                             error('Logic Error!');
%                         end
%                     else
%                         if single(dis(SL_location_coord(:,k),CL{5})+e)==single(e)
%                             if k==1
%                                 Zone_ID(1,k)=123;
%                                 A_matrix(k,:)=[0,0];
%                             else
%                                 error('Logic Error!');
%                             end
%                         else
%                             if single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor1)+e)==single(DIS_01+e)
%                                 Zone_ID(1,k)=31;
%                                 A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A2'*A2)*A2');
%                             elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor2)+e)==single(DIS_02+e)
%                                 Zone_ID(1,k)=12;
%                                 A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A3'*A3)*A3');
%                             elseif single(dis(SL_location_coord(:,k),CL{5})+dis(SL_location_coord(:,k),Coor3)+e)==single(DIS_03+e)
%                                 Zone_ID(1,k)=23;
%                                 A_matrix(k,:)=(SL_location_coord(:,k)-CL{5})'*(inv(A1'*A1)*A1');
%                             else
%                                 error('Logic Error!');
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             CL{64}=neigh_id;
%             CL{65}=[0,0,0];
%             CL{66}=A_matrix;
%             CL{67}=Zone_ID;
%             
%             CELL{i}=CL;
%         else
%             error('logic error!');
%         end
%     end
% end