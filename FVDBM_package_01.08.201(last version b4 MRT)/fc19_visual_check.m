d=((Y2-Y1)+(X2-X1))/10;
plot(X1-d,Y1-d,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X2+d,Y1-d,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X2+d,Y2+d,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
hold on
plot(X1-d,Y2+d,'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
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

for i=65:O
    FC=FACE{i};
    fc16=FC{16};
    fc17=FC{17};
    fc19=FC{19};
    fc23=FC{23};
    C_b=FC{7};
    if fc23==1
%     if FC{2}==1 % Check inner boundary
        %     if FC{2}==1 % Check outter boundary
        %     if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose further downwind stencil point is out of boundary
        %     if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose further upwind stencil point is out of boundary
        %     if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose both further downwind  and further upwind stencil points are out of boundary
        %     if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose all stencil points are within the boundary
        %% Plot the stencil
        p1=C_b-FC{4}'*d;
        p2=C_b+FC{4}'*d;
        plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
        pause(1)
        hold on
        %% Plot the further downwind stencil point
        if fc16(1,1)==0
            ;
        else
            c_dd=CELL{fc16(1,1)};
            C_dd=c_dd{5};
            plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
            pause(1)
            hold on;
        end
        plot(fc17(1,1),fc17(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
        pause(1)
        hold on;
        %% Plot the downwind stencil point
        if fc16(1,2)==0
            ;
        else
            c_d=CELL{fc16(1,2)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(1)
            hold on;
        end
        plot(fc17(1,2),fc17(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
        pause(1)
        hold on;
        %% Plot the upwind stencil point
        if fc16(1,3)==0
            ;
        else
            c_u=CELL{fc16(1,3)};
            C_u=c_u{5};
            plot(C_u(1),C_u(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(1)
            hold on;
        end
        plot(fc17(1,3),fc17(2,3),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
        pause(1)
        hold on;
        %% Plot the further upwind stencil point
        if fc16(1,4)==0
            ;
        else
            c_uu=CELL{fc16(1,4)};
            C_uu=c_uu{5};
            plot(C_uu(1),C_uu(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(1)
            hold on;
        end
        plot(fc17(1,4),fc17(2,4),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        pause(1)
        hold on;
        %% Plot the zone if the stencil point is located within the computational domain
        for j=1:length(fc16)
            if fc16(1,j)~=0
                cell_in=CELL{fc16(1,j)};
                C_nd1=cell_in{13};
                C_nd2=cell_in{14};
                C_nd3=cell_in{15};
                C_c=cell_in{5};
                if fc19(1,j)==1
                    edge1_x=[C_c(1,1),C_nd1(1,1)];
                    edge2_x=[C_nd1(1,1),C_nd2(1,1)];
                    edge3_x=[C_nd2(1,1),C_c(1,1)];
                    
                    edge1_y=[C_c(2,1),C_nd1(2,1)];
                    edge2_y=[C_nd1(2,1),C_nd2(2,1)];
                    edge3_y=[C_nd2(2,1),C_c(2,1)];
                elseif fc19(1,j)==2
                    edge1_x=[C_c(1,1),C_nd2(1,1)];
                    edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                    edge3_x=[C_nd3(1,1),C_c(1,1)];
                    
                    edge1_y=[C_c(2,1),C_nd2(2,1)];
                    edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                    edge3_y=[C_nd3(2,1),C_c(2,1)];
                elseif fc19(1,j)==3
                    edge1_x=[C_c(1,1),C_nd3(1,1)];
                    edge2_x=[C_nd3(1,1),C_nd1(1,1)];
                    edge3_x=[C_nd1(1,1),C_c(1,1)];
                    
                    edge1_y=[C_c(2,1),C_nd3(2,1)];
                    edge2_y=[C_nd3(2,1),C_nd1(2,1)];
                    edge3_y=[C_nd1(2,1),C_c(2,1)];
                elseif fc19(1,j)==4
                    edge1_x=[C_nd1(1,1),C_nd2(1,1)];
                    edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                    edge3_x=[C_nd3(1,1),C_nd1(1,1)];
                    
                    edge1_y=[C_nd1(2,1),C_nd2(2,1)];
                    edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                    edge3_y=[C_nd3(2,1),C_nd1(2,1)];
                else
                    error('Logic error!');
                end
                plot(edge1_x, edge1_y,'blue', 'linewidth',1);
                hold on;
                plot(edge2_x, edge2_y,'blue', 'linewidth',1);
                hold on;
                plot(edge3_x, edge3_y,'blue', 'linewidth',1);
                hold on;
                pause(1)
            end
        end
    end
end
