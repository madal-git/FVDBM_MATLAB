C_1=[4/3;1/3];
C_2=[10/3;2/3];
C_3=[5/3;11/3];
in_triangle(C_1,C_1,C_2,C_3)
Ori=[3;1/3];
A=[4/3;-1/3];
B=[2/3;1/3];
D=[2/3;2];
[~,C1]=intercept(Ori,(D-Ori),C_1,C_2);
[~,C2]=intercept(Ori,(D-Ori),C_1,C_3);
on_edge(C1,C_1,C_2)
on_edge(C2,C_1,C_3)

C=norm_joint(Ori,(B-Ori),A);
in_triangle(C,C_1,C_2,C_3)

in_triangle((C+Ori+b)/3,Ori,B,C)


on_edge(fc17(:,i),C_1,C_2)
on_edge(fc17(:,i),C_2,C_3)
on_edge(fc17(:,i),C_3,C_1)


    BL=[X1;Y1];
    BR=[X2;Y1];
    TL=[X1;Y2];
    TR=[X2;Y2];
    AA=(BL+BR)/2;
    on_edge(AA,BL,BR)
    
    
    C=CELL{1};
    C1=C{13};
    C2=C{14};
    C3=C{15};
    in_triangle(fc17(:,1),C1,C2,C3)
    on_edge(fc17(:,1),C1,C2)
    on_edge(fc17(:,1),C2,C3)
    on_edge(fc17(:,1),C1,C3)
    
    
    p=[1/3;1/2];
    T1=[1.23;0.4/3];
    T2=[-1/3;11/3];
    T3=[-5/3;-4/3];
    in_triangle(p,T1,T2,T3)
    edge1_x=[T1(1,1),T2(1,1)];
    edge2_x=[T2(1,1),T3(1,1)];
    edge3_x=[T3(1,1),T1(1,1)];
    
    edge1_y=[T1(2,1),T2(2,1)];
    edge2_y=[T2(2,1),T3(2,1)];
    edge3_y=[T3(2,1),T1(2,1)];
    
    plot(p(1),p(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
    hold on
    plot(edge1_x, edge1_y,'blue', 'linewidth',1);
    hold on;
    plot(edge2_x, edge2_y,'blue', 'linewidth',1);
    hold on;
    plot(edge3_x, edge3_y,'blue', 'linewidth',1);
    hold on;
    pause(1)
    
    % Dilation
    Factor=2;
    p_new=p*Factor;
    T1_new=T1*Factor;
    T2_new=T2*Factor;
    T3_new=T3*Factor;
    
    in_triangle(p_new,T1_new,T2_new,T3_new)
    edge1_x=[T1_new(1,1),T2_new(1,1)];
    edge2_x=[T2_new(1,1),T3_new(1,1)];
    edge3_x=[T3_new(1,1),T1_new(1,1)];
    
    edge1_y=[T1_new(2,1),T2_new(2,1)];
    edge2_y=[T2_new(2,1),T3_new(2,1)];
    edge3_y=[T3_new(2,1),T1_new(2,1)];
    
    plot(p_new(1),p_new(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
    hold on
    plot(edge1_x, edge1_y,'red', 'linewidth',1);
    hold on;
    plot(edge2_x, edge2_y,'red', 'linewidth',1);
    hold on;
    plot(edge3_x, edge3_y,'red', 'linewidth',1);
    hold on;