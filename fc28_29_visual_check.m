if FM==0
    T=0;
else
    T=0.5;
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
        fc25=FC{25};
        fc26=FC{26};
        fc27=FC{27};
        fc28=FC{28};
        fc29=FC{29};
        for i=1:length(fc24)
            Cell_fix=CELL{fc24(1,i)};
            C_fix=Cell_fix{5};
            %% Plot the stencil
            p1=C_b-FC{4}'*d;
            p2=C_b+FC{4}'*d;
            plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
            hold on
            if fc16(:,1)==0 && fc16(:,4)~=0
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
                pause(T)
                hold on
            elseif fc16(:,1)~=0 && fc16(:,4)==0
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
                pause(T)
                hold on
            elseif fc16(:,1)==0 && fc16(:,4)==0
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
                hold on
                
                
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
                pause(T)
                hold on
            else
                ;
            end
            % Plot stencil point and corresponding centroid
            
            c_dd=CELL{fc24(1,1)};
            C_dd=c_dd{5};
            plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
            pause(T)
            hold on;
            
            plot(fc25(1,1),fc25(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
            pause(T)
            hold on;
            
            c_d=CELL{fc24(1,2)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(T)
            hold on;
            
            plot(fc25(1,2),fc25(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
            pause(T)
            hold on;
            
            c_u=CELL{fc24(1,3)};
            C_u=c_u{5};
            plot(C_u(1),C_u(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(T)
            hold on;
            
            plot(fc25(1,3),fc25(2,3),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
            pause(T)
            hold on;
            
            c_uu=CELL{fc24(1,4)};
            C_uu=c_uu{5};
            plot(C_uu(1),C_uu(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(T)
            hold on;
            plot(fc25(1,4),fc25(2,4),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
            pause(T)
            hold on;
            break;
        end
        
        for i=1:length(fc24)
            Cell_fix=CELL{fc24(1,i)};
            C_fix=Cell_fix{5};
            
            C_ebcl=fc28(:,i);
            Cell1=CELL{C_ebcl(2,1)};
            Cell2=CELL{C_ebcl(3,1)};
            Cell3=CELL{C_ebcl(4,1)};
            C1=Cell1{5};
            C2=Cell2{5};
            C3=Cell3{5};
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
            C_nd1=fc29(1:2,i);
            C_nd2=fc29(3:4,i);
            C_nd3=fc29(5:6,i);
            edge1_x=[C_nd1(1,1),C_nd2(1,1)];
            edge2_x=[C_nd2(1,1),C_nd3(1,1)];
            edge3_x=[C_nd3(1,1),C_nd1(1,1)];
            
            edge1_y=[C_nd1(2,1),C_nd2(2,1)];
            edge2_y=[C_nd2(2,1),C_nd3(2,1)];
            edge3_y=[C_nd3(2,1),C_nd1(2,1)];
            
            plot(edge1_x, edge1_y,'red', 'linewidth',1);
            hold on;
            plot(edge2_x, edge2_y,'red', 'linewidth',1);
            hold on;
            plot(edge3_x, edge3_y,'red', 'linewidth',1);
            hold on;
            pause(T)
            
        end
        
        
    end
end
