plot(X1-0.5,Y1-0.5,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X2+0.5,Y1-0.5,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X2+0.5,Y2+0.5,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X1-0.5,Y2+0.5,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
axis equal tight;
for r=1:M;
    P=CELL{r};
    plot(P{22},P{23},'black', 'linewidth',2);
    hold on;
    plot(P{24},P{25},'black', 'linewidth',2);
    hold on;
    plot(P{26},P{27},'black', 'linewidth',2);
    hold on;
end;

% for i=1:O
%     FC=FACE{i};
%     fc18=FC{18};
%     fc19=FC{19};
%     fc20=FC{20};
%     C_b=FC{7};
%     
%     if fc19(1,1)==2 || fc19(1,2)==2
%         % Check
%         if FC{2}~=0
%             error('The face should be interior!');
%         end
%         %% Plot the face
%         p1=FC{5};
%         p2=FC{6};
%         plot([p1(1),p2(1)],[p1(2),p2(2)],'red', 'linewidth',2);
%         pause(0.1)
%         hold on
%         %% Plot the downwind stencil point
%         plot(fc18(1,1),fc18(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
%         pause(0.1)
%         hold on;
%         %% Plot the upwind stencil point
%         plot(fc18(1,2),fc18(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
%         pause(0.1)
%         hold on;
%         %% Plot the zone if the stencil point is located within the computational domain
%         for j=1:2
%                 if fc20(1,j)==0
%                     if fc20(2,j)==fc20(3,j) && fc20(3,j)==fc20(4,j)
%                         Cell1=CELL{fc20(2,j)};
%                         Cell2=CELL{fc20(3,j)};
%                         Cell3=CELL{fc20(4,j)};
%                         C_nd1=Cell1{13};
%                         C_nd2=Cell2{14};
%                         C_nd3=Cell3{15};
%                     else
%                         Cell1=CELL{fc20(2,j)};
%                         Cell2=CELL{fc20(3,j)};
%                         Cell3=CELL{fc20(4,j)};
%                         C_nd1=Cell1{5};
%                         C_nd2=Cell2{5};
%                         C_nd3=Cell3{5};
%                     end
%                 elseif fc20(1,j)==1
%                     Node1=NODE{fc20(2,j)};
%                     Cell2=CELL{fc20(3,j)};
%                     Cell3=CELL{fc20(4,j)};
%                     C_nd1=Node1{3};
%                     C_nd2=Cell2{5};
%                     C_nd3=Cell3{5};
%                 elseif fc20(1,j)==2
%                     Node1=NODE{fc20(2,j)};
%                     Node2=NODE{fc20(3,j)};
%                     Cell3=CELL{fc20(4,j)};
%                     C_nd1=Node1{3};
%                     C_nd2=Node2{3};
%                     C_nd3=Cell3{5};
%                 else
%                     error('The points that enlose the stencil point cannot be all nodes!');
%                 end
%                 
%                 edge1_x=[C_nd1(1,1),C_nd2(1,1)];
%                 edge2_x=[C_nd2(1,1),C_nd3(1,1)];
%                 edge3_x=[C_nd3(1,1),C_nd1(1,1)];
%                 
%                 edge1_y=[C_nd1(2,1),C_nd2(2,1)];
%                 edge2_y=[C_nd2(2,1),C_nd3(2,1)];
%                 edge3_y=[C_nd3(2,1),C_nd1(2,1)];
%                 
%                 plot(edge1_x, edge1_y,'blue', 'linewidth',1);
%                 hold on;
%                 plot(edge2_x, edge2_y,'blue', 'linewidth',1);
%                 hold on;
%                 plot(edge3_x, edge3_y,'blue', 'linewidth',1);
%                 hold on;
%                 pause(0.1)
%         end
%     end
% end

% for i=1:O
%     FC=FACE{i};
%     if FC{2}==1
%         neigh_up_p=FC{14};
%         neigh_down_p=FC{15};
%         Cell_up_p=CELL{neigh_up_p(1)};
%         Cell_down_p=CELL{neigh_down_p(1)};
%         C_u=Cell_up_p{5};
%         C_d=Cell_down_p{5};
%         plot(C_u(1,1),C_u(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
%         hold on;
%         plot(C_d(1,1),C_d(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
%         hold on;
%     end
% end

for i=1:N
    ND=NODE{i};
    if ND{2}==1
        L=ND{8};
        Pool=ND{9};
        for j=1:L
            Cell_s=CELL{Pool(j)};
            C_s=Cell_s{5};
            plot(C_s(1,1),C_s(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
            hold on;
        end
        pause(1)
    end
end
