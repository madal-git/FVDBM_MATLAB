if FM==0
    T=0.6;
else
    T=1;
end
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

for l=1:O
    FC=FACE{l};
    fc16=FC{16};
    if FC{23}==4 %&& length(union(0,fc16))<length(fc16)
        C_b=FC{7};
        fc18=FC{18};
        fc22=FC{22};
        fc24=FC{24};
        fc26=FC{26};
        fc27=FC{27};
        fc28=FC{28};
        fc29=FC{29};
        
        fc37=FC{37};
        fc38=FC{38};
        fc40=FC{40};
        fc42=FC{42};
        fc43=FC{43};
        p1=C_b-FC{4}'*d;
        p2=C_b+FC{4}'*d;
        plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
        hold on
        if fc16(:,1)==0 && fc16(:,4)~=0
            if fc37(1,1)~=0
                Nd1=NODE{fc18(1,1)};
                Nd2=NODE{fc18(2,1)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc18(3,1)*L*v;
                if fc22(1,1)==1
                    C_i_m=C_i+[0;-(Y2-Y1)];
                elseif fc22(1,1)==2
                    C_i_m=C_i+[-(X2-X1);0];
                elseif fc22(1,1)==3
                    C_i_m=C_i+[0;(Y2-Y1)];
                elseif fc22(1,1)==4
                    C_i_m=C_i+[(X2-X1);0];
                else
                    error('Logic error!');
                end
                p1=C_i_m-FC{4}'*d;
                p2=C_i_m+FC{4}'*d;
                plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
                pause(1)
                hold on
            else
                pause(1)
            end
        elseif fc16(:,1)~=0 && fc16(:,4)==0
            if fc37(1,4)~=0
                Nd1=NODE{fc18(1,4)};
                Nd2=NODE{fc18(2,4)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc18(3,4)*L*v;
                if fc22(1,4)==1
                    C_i_m=C_i+[0;-(Y2-Y1)];
                elseif fc22(1,4)==2
                    C_i_m=C_i+[-(X2-X1);0];
                elseif fc22(1,4)==3
                    C_i_m=C_i+[0;(Y2-Y1)];
                elseif fc22(1,4)==4
                    C_i_m=C_i+[(X2-X1);0];
                else
                    error('Logic error!');
                end
                p1=C_i_m-FC{4}'*d;
                p2=C_i_m+FC{4}'*d;
                plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
                pause(1)
                hold on
            else
                pause(1)
            end
        elseif fc16(:,1)==0 && fc16(:,4)==0
            if fc37(1,1)~=0
                Nd1_1=NODE{fc18(1,1)};
                Nd2_1=NODE{fc18(2,1)};
                L=dis(Nd1_1{3},Nd2_1{3});
                v=Nd2_1{3}-Nd1_1{3};
                v=v/norm(v);
                C_i_1=Nd1_1{3}+fc18(3,1)*L*v;
                if fc22(1,1)==1
                    C_i_m_1=C_i_1+[0;-(Y2-Y1)];
                elseif fc22(1,1)==2
                    C_i_m_1=C_i_1+[-(X2-X1);0];
                elseif fc22(1,1)==3
                    C_i_m_1=C_i_1+[0;(Y2-Y1)];
                elseif fc22(1,1)==4
                    C_i_m_1=C_i_1+[(X2-X1);0];
                else
                    error('Logic error!');
                end
                p1=C_i_m_1-FC{4}'*d;
                p2=C_i_m_1+FC{4}'*d;
                plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
                pause(1)
                hold on
            else
                pause(1)
            end
            
            if fc37(1,4)~=0
                Nd1_2=NODE{fc18(1,4)};
                Nd2_2=NODE{fc18(2,4)};
                L=dis(Nd1_2{3},Nd2_2{3});
                v=Nd2_2{3}-Nd1_2{3};
                v=v/norm(v);
                C_i_2=Nd1_2{3}+fc18(3,4)*L*v;
                if fc22(1,4)==1
                    C_i_m_2=C_i_2+[0;-(Y2-Y1)];
                elseif fc22(1,4)==2
                    C_i_m_2=C_i_2+[-(X2-X1);0];
                elseif fc22(1,4)==3
                    C_i_m_2=C_i_2+[0;(Y2-Y1)];
                elseif fc22(1,4)==4
                    C_i_m_2=C_i_2+[(X2-X1);0];
                else
                    error('Logic error!');
                end
                p1=C_i_m_2-FC{4}'*d;
                p2=C_i_m_2+FC{4}'*d;
                plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
                pause(1)
                hold on
            else
                pause(1)
            end
        else
            ;
        end
        %% Plot the further downwind stencil point
        if fc37(1,1)==0
            % Plot the intercept point on boundary
            Nd1=NODE{fc40(1,1)};
            Nd2=NODE{fc40(2,1)};
            L=dis(Nd1{3},Nd2{3});
            v=Nd2{3}-Nd1{3};
            v=v/norm(v);
            C_i=Nd1{3}+fc40(3,1)*L*v;
            plot(C_i(1),C_i(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
            pause(1)
            hold on;
        else
            c_dd=CELL{fc37(1,1)};
            C_dd=c_dd{5};
            plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
            pause(1)
            hold on;
        end
        plot(fc38(1,1),fc38(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
        pause(1)
        hold on;
        %% Plot the downwind stencil point
        if fc37(1,2)==0
            if fc37(1,1)==0
                ;
            else
                % Plot the intercept point on boundary
                Nd1=NODE{fc40(1,2)};
                Nd2=NODE{fc40(2,2)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc40(3,2)*L*v;
                plot(C_i(1),C_i(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
                pause(1)
                hold on;
            end
        else
            c_d=CELL{fc37(1,2)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(1)
            hold on;
        end
        plot(fc38(1,2),fc38(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
        pause(1)
        hold on;
        %% Plot the upwind stencil point
        if fc37(1,3)==0
            % Plot the intercept point on boundary
            Nd1=NODE{fc40(1,3)};
            Nd2=NODE{fc40(2,3)};
            L=dis(Nd1{3},Nd2{3});
            v=Nd2{3}-Nd1{3};
            v=v/norm(v);
            C_i=Nd1{3}+fc40(3,3)*L*v;
            plot(C_i(1),C_i(2),'Marker', 'x','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
            pause(1)
            hold on;
        else
            c_u=CELL{fc37(1,3)};
            C_u=c_u{5};
            plot(C_u(1),C_u(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(1)
            hold on;
        end
        plot(fc38(1,3),fc38(2,3),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
        pause(1)
        hold on;
        %% Plot the further upwind stencil point
        if fc37(1,4)==0
            % Plot the intercept point on boundary
            if fc37(1,3)==0
                ;
            else
                Nd1=NODE{fc40(1,4)};
                Nd2=NODE{fc40(2,4)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc40(3,4)*L*v;
                plot(C_i(1),C_i(2),'Marker', 'x','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
                pause(1)
                hold on;
            end
        else
            c_uu=CELL{fc37(1,4)};
            C_uu=c_uu{5};
            plot(C_uu(1),C_uu(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(1)
            hold on;
        end
        plot(fc38(1,4),fc38(2,4),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        pause(1)
        hold on;
        
        %% Plot enclosing points
        for i=1:length(fc37)
            if fc37(1,i)~=0
                C_ebcl=fc42(:,i);
                if fc42(1,i)==0
                    Cell1=CELL{C_ebcl(2,1)};
                    Cell2=CELL{C_ebcl(3,1)};
                    Cell3=CELL{C_ebcl(4,1)};
                    if fc42(2,i)==fc42(3,i) && fc42(3,i)==fc42(4,i)
                        C1=Cell1{13};
                        C2=Cell2{14};
                        C3=Cell3{15};
                    else
                        C1=Cell1{5};
                        C2=Cell2{5};
                        C3=Cell3{5};
                    end
                elseif fc42(1,i)==1
                    Node1=NODE{C_ebcl(2,1)};
                    Cell2=CELL{C_ebcl(3,1)};
                    Cell3=CELL{C_ebcl(4,1)};
                    C1=Node1{3};
                    C2=Cell2{5};
                    C3=Cell3{5};
                elseif fc42(1,i)==2
                    Node1=NODE{C_ebcl(2,1)};
                    Node2=NODE{C_ebcl(3,1)};
                    Cell3=CELL{C_ebcl(4,1)};
                    C1=Node1{3};
                    C2=Node2{3};
                    C3=Cell3{5};
                else
                    error('adf');
                end
                
                if i==1
                    plot(C1(1),C1(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'c');
                    pause(T)
                    hold on;
                    plot(C2(1),C2(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'c');
                    pause(T)
                    hold on;
                    plot(C3(1),C3(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'c');
                    pause(T)
                    hold on;
                elseif i==2
                    plot(C1(1),C1(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'g');
                    pause(T)
                    hold on;
                    plot(C2(1),C2(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'g');
                    pause(T)
                    hold on;
                    plot(C3(1),C3(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'g');
                    pause(T)
                    hold on;
                elseif i==3
                    plot(C1(1),C1(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'm');
                    pause(T)
                    hold on;
                    plot(C2(1),C2(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'm');
                    pause(T)
                    hold on;
                    plot(C3(1),C3(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'm');
                    pause(T)
                    hold on;
                elseif i==4
                    plot(C1(1),C1(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'r');
                    pause(T)
                    hold on;
                    plot(C2(1),C2(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'r');
                    pause(T)
                    hold on;
                    plot(C3(1),C3(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'r');
                    pause(T)
                    hold on;
                else
                    ;
                end
                C_nd1=fc43(1:2,i);
                C_nd2=fc43(3:4,i);
                C_nd3=fc43(5:6,i);
                edge1_x=[C_nd1(1,1),C_nd2(1,1)];
                edge2_x=[C_nd2(1,1),C_nd3(1,1)];
                edge3_x=[C_nd3(1,1),C_nd1(1,1)];
                
                edge1_y=[C_nd1(2,1),C_nd2(2,1)];
                edge2_y=[C_nd2(2,1),C_nd3(2,1)];
                edge3_y=[C_nd3(2,1),C_nd1(2,1)];
                if i==1
                    plot(edge1_x, edge1_y,'c', 'linewidth',1);
                    hold on;
                    plot(edge2_x, edge2_y,'c', 'linewidth',1);
                    hold on;
                    plot(edge3_x, edge3_y,'c', 'linewidth',1);
                    hold on;
                    pause(T)
                elseif i==2
                    plot(edge1_x, edge1_y,'g', 'linewidth',1);
                    hold on;
                    plot(edge2_x, edge2_y,'g', 'linewidth',1);
                    hold on;
                    plot(edge3_x, edge3_y,'g', 'linewidth',1);
                    hold on;
                    pause(T)
                elseif i==3
                    plot(edge1_x, edge1_y,'m', 'linewidth',1);
                    hold on;
                    plot(edge2_x, edge2_y,'m', 'linewidth',1);
                    hold on;
                    plot(edge3_x, edge3_y,'m', 'linewidth',1);
                    hold on;
                    pause(T)
                elseif i==4
                    plot(edge1_x, edge1_y,'red', 'linewidth',1);
                    hold on;
                    plot(edge2_x, edge2_y,'red', 'linewidth',1);
                    hold on;
                    plot(edge3_x, edge3_y,'red', 'linewidth',1);
                    hold on;
                    pause(T)
                else
                    ;
                end
                
            end
        end
        
        
    end
end
