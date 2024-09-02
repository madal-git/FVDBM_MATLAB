if FM==0
    T=0.6;
else
    T=0.6;
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

ND=cell(1,4);
ND{1}=NODE{N_I+N_L-1};
ND{2}=NODE{N_I+N_L-1+N_H-1};
ND{3}=NODE{N_I+N_L-1+N_H-1+N_L-1};
ND{4}=NODE{N};

for l=1:4
    Nd=ND{l};
    C=Nd{3};
    if Nd{2}<73
        Pool=Nd{8};
        Nc=Nd{9};
    else
        Pool=Nd{4};
        Nc=Nd{5};
    end
    for i=1:Pool
        if l==1
            plot(C(1),C(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'k');
            pause(T)
            hold on;
        elseif l==2
            plot(C(1),C(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'b','MarkerEdgeColor', 'k');
            pause(T)
            hold on;
        elseif l==3
            plot(C(1),C(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'k');
            pause(T)
            hold on;
        else
            plot(C(1),C(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'k');
            pause(T)
            hold on;
        end
        P=CELL{Nc(i)};
        C_nd1=P{13};
        C_nd2=P{14};
        C_nd3=P{15};
        edge1_x=[C_nd1(1,1),C_nd2(1,1)];
        edge2_x=[C_nd2(1,1),C_nd3(1,1)];
        edge3_x=[C_nd3(1,1),C_nd1(1,1)];
        
        edge1_y=[C_nd1(2,1),C_nd2(2,1)];
        edge2_y=[C_nd2(2,1),C_nd3(2,1)];
        edge3_y=[C_nd3(2,1),C_nd1(2,1)];
        if l==1
            plot(edge1_x, edge1_y,'red', 'linewidth',1);
            hold on;
            plot(edge2_x, edge2_y,'red', 'linewidth',1);
            hold on;
            plot(edge3_x, edge3_y,'red', 'linewidth',1);
            hold on;
            pause(T)
        elseif l==2
            plot(edge1_x, edge1_y,'b', 'linewidth',1);
            hold on;
            plot(edge2_x, edge2_y,'b', 'linewidth',1);
            hold on;
            plot(edge3_x, edge3_y,'b', 'linewidth',1);
            hold on;
            pause(T)
        elseif l==3
            plot(edge1_x, edge1_y,'m', 'linewidth',1);
            hold on;
            plot(edge2_x, edge2_y,'m', 'linewidth',1);
            hold on;
            plot(edge3_x, edge3_y,'m', 'linewidth',1);
            hold on;
            pause(T)
        else
            plot(edge1_x, edge1_y,'c', 'linewidth',1);
            hold on;
            plot(edge2_x, edge2_y,'c', 'linewidth',1);
            hold on;
            plot(edge3_x, edge3_y,'c', 'linewidth',1);
            hold on;
            pause(T)
        end
    end
end