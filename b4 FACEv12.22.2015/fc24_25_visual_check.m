if FM==0
    T=1;
else
    T=2;
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
    if FC{2}>=0 % Exclude the inner boundary
        fc16=FC{16};
        fc22=FC{22};
%         if length(union(0,fc16))<=length(fc16) && (fc22(1,1)>=0 && fc22(1,end)>=0)% Face type 2 and 3, which has at least one of the stencil point is out of boundary
        if length(union(0,fc16))==length(fc16) && (fc22(1,1)>=0 && fc22(1,end)>=0)% Face type 2 and 3, which has at least one of the stencil point is out of boundary
            C_b=FC{7};
            fc18=FC{18};
            fc24=FC{24};
            fc25=FC{25};
            fc26=FC{26};
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
                error('Logic error!');
            end

            
            %% Plot the further downwind stencil point
            c_dd=CELL{fc24(1,1)};
            C_dd=c_dd{5};
            plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
            pause(T)
            hold on;

            plot(fc25(1,1),fc25(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
            pause(T)
            hold on;
            %% Plot the further downwind stencil point
            c_d=CELL{fc24(1,2)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(T)
            hold on;
            
            plot(fc25(1,2),fc25(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
            pause(T)
            hold on;
            %% Plot the upwind stencil point
            c_u=CELL{fc24(1,3)};
            C_u=c_u{5};
            plot(C_u(1),C_u(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(T)
            hold on;
            
            plot(fc25(1,3),fc25(2,3),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
            pause(T)
            hold on;
            %% Plot the further upwind stencil point
            c_uu=CELL{fc24(1,4)};
            C_uu=c_uu{5};
            plot(C_uu(1),C_uu(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(T)
            hold on;
            
            plot(fc25(1,4),fc25(2,4),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
            pause(T)
            hold on;
        end
    end
end
