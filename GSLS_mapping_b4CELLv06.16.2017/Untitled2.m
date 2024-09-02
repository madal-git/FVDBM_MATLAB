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
l=201;
FC=FACE{l};
fc16=FC{16};
fc17=FC{17};
fc18=FC{18};
fc23=FC{23};
C_b=FC{7};

i=4;
C=CELL{fc16(1,i)};
C_fix=C{5};
C_nd1=C{13};
C_nd2=C{14};
C_nd3=C{15};

p1=C_b-FC{4}'*d;
p2=C_b+FC{4}'*d;
plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
pause(1)
hold on

c_dd=CELL{fc16(1,i)};
C_dd=c_dd{5};
plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
pause(1)
hold on;

plot(fc17(1,i),fc17(2,i),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
pause(1)
hold on;

edge1_x=[C_nd1(1,1),C_nd2(1,1)];
edge2_x=[C_nd2(1,1),C_nd3(1,1)];
edge3_x=[C_nd3(1,1),C_nd1(1,1)];

edge1_y=[C_nd1(2,1),C_nd2(2,1)];
edge2_y=[C_nd2(2,1),C_nd3(2,1)];
edge3_y=[C_nd3(2,1),C_nd1(2,1)];

plot(edge1_x, edge1_y,'blue', 'linewidth',1);
hold on;
plot(edge2_x, edge2_y,'blue', 'linewidth',1);
hold on;
plot(edge3_x, edge3_y,'blue', 'linewidth',1);
hold on;
pause(2)

edg1_x=[C_fix(1,1),C_nd1(1,1)];
edg2_x=[C_fix(1,1),C_nd2(1,1)];
edg3_x=[C_fix(1,1),C_nd3(1,1)];

edg1_y=[C_fix(2,1),C_nd1(2,1)];
edg2_y=[C_fix(2,1),C_nd2(2,1)];
edg3_y=[C_fix(2,1),C_nd3(2,1)];

plot(edg1_x, edg1_y,'red', 'linewidth',1);
hold on;
plot(edg2_x, edg2_y,'red', 'linewidth',1);
hold on;
plot(edg3_x, edg3_y,'red', 'linewidth',1);
hold on;
pause(2)

in_triangle(fc17(:,i),C_fix,C_nd1,C_nd2)
in_triangle(fc17(:,i),C_fix,C_nd2,C_nd3)
in_triangle(fc17(:,i),C_fix,C_nd3,C_nd1)

ND=NODE{80252};
C_nd=ND{3};
C1=CELL{207};
C_1=C1{5};
C2=CELL{1004};
C_2=C2{5};
Cnd1=C2{13};
Cnd2=C2{14};
Cnd3=C2{15};
ed1_x=[C_nd(1,1),C_1(1,1)];
ed2_x=[C_1(1,1),C_2(1,1)];
ed3_x=[C_2(1,1),C_nd(1,1)];

ed1_y=[C_nd(2,1),C_1(2,1)];
ed2_y=[C_1(2,1),C_2(2,1)];
ed3_y=[C_2(2,1),C_nd(2,1)];
plot(ed1_x, ed1_y,'k', 'linewidth',1);
hold on;
plot(ed2_x, ed2_y,'k', 'linewidth',1);
hold on;
plot(ed3_x, ed3_y,'k', 'linewidth',1);
hold on;
pause(2)

e1_x=[Cnd1(1,1),Cnd2(1,1)];
e2_x=[Cnd2(1,1),Cnd3(1,1)];
e3_x=[Cnd3(1,1),Cnd1(1,1)];

e1_y=[Cnd1(2,1),Cnd2(2,1)];
e2_y=[Cnd2(2,1),Cnd3(2,1)];
e3_y=[Cnd3(2,1),Cnd1(2,1)];
plot(e1_x, e1_y,'k', 'linewidth',1);
hold on;
plot(e2_x, e2_y,'k', 'linewidth',1);
hold on;
plot(e3_x, e3_y,'k', 'linewidth',1);
hold on;
pause(2)

Cell_fix=C;
C_fixed=C{5};
fc19=FC{19};
Face_target=FACE{Cell_fix{15+fc19(1,i)}};
if Face_target{10}~=0 && Face_target{11}==0
    ND=NODE{Face_target{8}};
    C_nd=ND{3};
    nd=Face_target{8};
elseif Face_target{10}==0 && Face_target{11}~=0
    ND=NODE{Face_target{8}};
    C_nd=ND{3};
    nd=Face_target{9};
else
    error('There must be one node on the boundary and the other is not!');
end
ND1=NODE{Face_target{8}};
ND2=NODE{Face_target{9}};
Cell_third_pool=union(ND1{5},ND2{5});
Cell_third_pool=setxor(fc16(1,i),Cell_third_pool);
L=length(Cell_third_pool);
a=0;
in_tri=0;
for k=1:L
    Cell_third=CELL{Cell_third_pool(k)};
    C_third_cell=Cell_third{5};
    if in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
        a=a+1;
        in_tri(a)=k;
    end
end

for j=1:a
    P=CELL{Cell_third_pool(in_tri(j))};
    C_C=P{5};
    CND1=P{13};
    CND2=P{14};
    CND3=P{15};
    plot(C_C(1),C_C(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
    pause(1)
    hold on;
    e1_x=[CND1(1,1),CND2(1,1)];
    e2_x=[CND2(1,1),CND3(1,1)];
    e3_x=[CND3(1,1),CND1(1,1)];
    
    e1_y=[CND1(2,1),CND2(2,1)];
    e2_y=[CND2(2,1),CND3(2,1)];
    e3_y=[CND3(2,1),CND1(2,1)];
    plot(e1_x, e1_y,'r', 'linewidth',1);
    hold on;
    plot(e2_x, e2_y,'r', 'linewidth',1);
    hold on;
    plot(e3_x, e3_y,'r', 'linewidth',1);
    hold on;
    pause(2)
end

fc20=FC{20};
Node1=NODE{fc20(2,i)};
Cell2=CELL{fc20(3,i)};
Cell3=CELL{fc20(4,i)};
in_triangle(fc17(:,i),Node1{3},Cell2{5},Cell3{5})
% on_edge()
if a==0
    error('Enclosing triangle is not found!');
end
%% Find the cell that has the shortest distance to the stencil point
Dis=zeros(1,a);
for k=1:a
    Cell_third=CELL{Cell_third_pool(in_tri(k))};
    C_third_cell=Cell_third{5};
    Dis(1,k)=dis(fc17(:,i),C_third_cell);
end
Dis_min=min(Dis);
for k=1:a
    if single(e+Dis_min)==single(e+Dis(1,k))
        break;
    end
end
if k==a
    if Dis(1,k)~=Dis_min
        error('Logic error!');
    end
end
% Third_cell_found=Cell_third_pool(in_tri(k));
% % % Check
Cell_third=CELL{Cell_third_pool(7)};
C_third_cell=Cell_third{5};
in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
on_edge(fc17(:,i),C_fixed,C_nd)
on_edge(fc17(:,i),C_fixed,C_third_cell)
on_edge(fc17(:,i),C_third_cell,C_nd)