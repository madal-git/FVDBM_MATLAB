VV=V2;

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

for i=1:O
    FC=FACE{i};
    fc16=FC{16};
    fc17=FC{17};
    fc18=FC{18};
    fc23=FC{23};
    S=FC{32};
    C_b=FC{7};
    if fc23==1
    %     if FC{2}==-1 % Check inner boundary
%         if FC{2}==1 % Check outter boundary
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
        s1=S{1};
        s2=S{2};
        s3=S{3};
        s4=S{4};
        s5=S{5};
        s6=S{6};
        s7=S{7};
        s8=S{8};
        s9=S{9};
        L=length(s1);
        for j=1:L
            %% Plot lattice velocity
            quiver(C_b(1,1),C_b(2,1),VV(1,j),VV(2,j),d/2);
            pause(1)
            hold on
            %% Plot the downwind centroid and stencil point, and circling points
            c_d=CELL{s1(1,j)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(1)
            hold on;
            
            plot(s4(1,j),s4(2,j),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
            pause(1)
            hold on;
            if s1(1,j)~=0
                if s7(1,j)==0
                    if s7(2,j)==s7(3,j) && s7(3,j)==s7(4,j)
                        Cell1=CELL{s7(2,j)};
                        Cell2=CELL{s7(3,j)};
                        Cell3=CELL{s7(4,j)};
                        C_nd1=Cell1{13};
                        C_nd2=Cell2{14};
                        C_nd3=Cell3{15};
                    else
                        Cell1=CELL{s7(2,j)};
                        Cell2=CELL{s7(3,j)};
                        Cell3=CELL{s7(4,j)};
                        C_nd1=Cell1{5};
                        C_nd2=Cell2{5};
                        C_nd3=Cell3{5};
                    end
                elseif s7(1,j)==1
                    Node1=NODE{s7(2,j)};
                    Cell2=CELL{s7(3,j)};
                    Cell3=CELL{s7(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Cell2{5};
                    C_nd3=Cell3{5};
                elseif s7(1,j)==2
                    Node1=NODE{s7(2,j)};
                    Node2=NODE{s7(3,j)};
                    Cell3=CELL{s7(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Node2{3};
                    C_nd3=Cell3{5};
                else
                    error('The stencil point cannot be enclosed by all boundary nodes!');
                end
                
                edge1_x=[C_nd1(1,1),C_nd2(1,1)];
                edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                edge3_x=[C_nd3(1,1),C_nd1(1,1)];
                
                edge1_y=[C_nd1(2,1),C_nd2(2,1)];
                edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                edge3_y=[C_nd3(2,1),C_nd1(2,1)];
                
                plot(edge1_x, edge1_y,'g', 'linewidth',1);
                hold on;
                plot(edge2_x, edge2_y,'g', 'linewidth',1);
                hold on;
                plot(edge3_x, edge3_y,'g', 'linewidth',1);
                hold on;
                pause(1)
            end
            %% Plot the upwind centroid and stencil point, and circling points
            c_d=CELL{s2(1,j)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(1)
            hold on;
            
            plot(s5(1,j),s5(2,j),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
            pause(1)
            hold on;
            
            if s2(1,j)~=0
                if s8(1,j)==0
                    if s8(2,j)==s8(3,j) && s8(3,j)==s8(4,j)
                        Cell1=CELL{s8(2,j)};
                        Cell2=CELL{s8(3,j)};
                        Cell3=CELL{s8(4,j)};
                        C_nd1=Cell1{13};
                        C_nd2=Cell2{14};
                        C_nd3=Cell3{15};
                    else
                        Cell1=CELL{s8(2,j)};
                        Cell2=CELL{s8(3,j)};
                        Cell3=CELL{s8(4,j)};
                        C_nd1=Cell1{5};
                        C_nd2=Cell2{5};
                        C_nd3=Cell3{5};
                    end
                elseif s8(1,j)==1
                    Node1=NODE{s8(2,j)};
                    Cell2=CELL{s8(3,j)};
                    Cell3=CELL{s8(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Cell2{5};
                    C_nd3=Cell3{5};
                elseif s8(1,j)==2
                    Node1=NODE{s8(2,j)};
                    Node2=NODE{s8(3,j)};
                    Cell3=CELL{s8(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Node2{3};
                    C_nd3=Cell3{5};
                else
                    error('The stencil point cannot be enclosed by all boundary nodes!');
                end
                
                edge1_x=[C_nd1(1,1),C_nd2(1,1)];
                edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                edge3_x=[C_nd3(1,1),C_nd1(1,1)];
                
                edge1_y=[C_nd1(2,1),C_nd2(2,1)];
                edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                edge3_y=[C_nd3(2,1),C_nd1(2,1)];
                
                plot(edge1_x, edge1_y,'m', 'linewidth',1);
                hold on;
                plot(edge2_x, edge2_y,'m', 'linewidth',1);
                hold on;
                plot(edge3_x, edge3_y,'m', 'linewidth',1);
                hold on;
                pause(1)
            end
            %% Plot the further upwind centroid and stencil point, and circling points
            c_d=CELL{s3(1,j)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(1)
            hold on;
            
            plot(s6(1,j),s6(2,j),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
            pause(1)
            hold on;
            
            if s3(1,j)~=0
                if s9(1,j)==0
                    if s9(2,j)==s9(3,j) && s9(3,j)==s9(4,j)
                        Cell1=CELL{s9(2,j)};
                        Cell2=CELL{s9(3,j)};
                        Cell3=CELL{s9(4,j)};
                        C_nd1=Cell1{13};
                        C_nd2=Cell2{14};
                        C_nd3=Cell3{15};
                    else
                        Cell1=CELL{s9(2,j)};
                        Cell2=CELL{s9(3,j)};
                        Cell3=CELL{s9(4,j)};
                        C_nd1=Cell1{5};
                        C_nd2=Cell2{5};
                        C_nd3=Cell3{5};
                    end
                elseif s9(1,j)==1
                    Node1=NODE{s9(2,j)};
                    Cell2=CELL{s9(3,j)};
                    Cell3=CELL{s9(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Cell2{5};
                    C_nd3=Cell3{5};
                elseif s9(1,j)==2
                    Node1=NODE{s9(2,j)};
                    Node2=NODE{s9(3,j)};
                    Cell3=CELL{s9(4,j)};
                    C_nd1=Node1{3};
                    C_nd2=Node2{3};
                    C_nd3=Cell3{5};
                else
                    error('The stencil point cannot be enclosed by all boundary nodes!');
                end
                
                edge1_x=[C_nd1(1,1),C_nd2(1,1)];
                edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                edge3_x=[C_nd3(1,1),C_nd1(1,1)];
                
                edge1_y=[C_nd1(2,1),C_nd2(2,1)];
                edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                edge3_y=[C_nd3(2,1),C_nd1(2,1)];
                
                plot(edge1_x, edge1_y,'r', 'linewidth',1);
                hold on;
                plot(edge2_x, edge2_y,'r', 'linewidth',1);
                hold on;
                plot(edge3_x, edge3_y,'r', 'linewidth',1);
                hold on;
                pause(1)
            end
        end
    end
end
