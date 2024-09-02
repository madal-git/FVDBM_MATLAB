plt_mesh;
hold on

for l=1:N
    ND=NODE{l};
    if ND{2}~=0
        nd22=ND{22};
        nd23=ND{23};
        nd24=ND{24};
        nd25=ND{25};
        
        %% plot the boundary node
        nd_bc_coord=ND{3};
        plot(nd_bc_coord(1),nd_bc_coord(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        pause(0.5)
        hold on
        
        %% plot the boundary normal vector
        n_end_coord=nd_bc_coord+nd22;
        n_x=[nd_bc_coord(1,1),n_end_coord(1,1)];
        n_y=[nd_bc_coord(2,1),n_end_coord(2,1)];
        plot(n_x, n_y,'r', 'linewidth',1);
        pause(0.5)
        hold on
        
        %% plot the centroid of the cell that contains the nearer stencil point
        CL_near=CELL{nd23(1,1)};
        centroid_cell_near=CL_near{5};
        plot(centroid_cell_near(1),centroid_cell_near(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'blue','MarkerEdgeColor', 'blue');
        pause(0.5)
        hold on
        
        %% plot the nearer stencil point
        plot(nd24(1,1),nd24(2,1),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'blue');
        pause(0.5)
        hold on
        
        %% plot the centroid of the cell that contains the further stencil point
        CL_further=CELL{nd23(1,2)};
        centroid_cell_further=CL_further{5};
        plot(centroid_cell_further(1),centroid_cell_further(2),'Marker', 's','Markersize',4, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
        pause(0.5)
        hold on
        
        %% plot the further stencil point
        plot(nd24(1,2),nd24(2,2),'Marker', 's','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'black');
        pause(0.5)
        hold on
        
    end
end