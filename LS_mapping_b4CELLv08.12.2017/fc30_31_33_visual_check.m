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
    fc18=FC{18};
    fc22=FC{22};
    fc30=FC{30};
    fc31=FC{31};
    fc33=FC{33};
    fc23=FC{23};
    C_b=FC{7};
    if fc23==2
        %     if FC{2}==-1 % Check inner boundary
        %         if FC{2}==1 % Check outter boundary
        %     if (fc30(1,1)==0 && fc30(1,2)~=0) && (fc30(1,3)~=0 && fc30(1,4)~=0) % Check the face whose further downwind stencil point is out of boundary
        %     if (fc30(1,1)~=0 && fc30(1,2)~=0) && (fc30(1,3)~=0 && fc30(1,4)==0) % Check the face whose further upwind stencil point is out of boundary
        %     if (fc30(1,1)==0 && fc30(1,2)~=0) && (fc30(1,3)~=0 && fc30(1,4)==0) % Check the face whose both further downwind  and further upwind stencil points are out of boundary
        %     if (fc30(1,1)~=0 && fc30(1,2)~=0) && (fc30(1,3)~=0 && fc30(1,4)~=0) % Check the face whose all stencil points are within the boundary
        %% Plot the stencil
        p1=C_b-FC{4}'*d;
        p2=C_b+FC{4}'*d;
        plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
        hold on
        if fc16(:,1)==0 && fc16(:,4)~=0
            if fc30(1,1)~=0
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
            if fc30(1,4)~=0
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
            if fc30(1,1)~=0
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
            
            if fc30(1,4)~=0
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
        if fc30(1,1)==0
            % Plot the intercept point on boundary
            Nd1=NODE{fc33(1,1)};
            Nd2=NODE{fc33(2,1)};
            L=dis(Nd1{3},Nd2{3});
            v=Nd2{3}-Nd1{3};
            v=v/norm(v);
            C_i=Nd1{3}+fc33(3,1)*L*v;
            plot(C_i(1),C_i(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
            pause(1)
            hold on;
        else
            c_dd=CELL{fc30(1,1)};
            C_dd=c_dd{5};
            plot(C_dd(1),C_dd(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'c');
            pause(1)
            hold on;
        end
        plot(fc31(1,1),fc31(2,1),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'c','MarkerEdgeColor', 'c');
        pause(1)
        hold on;
        %% Plot the downwind stencil point
        if fc30(1,2)==0
            if fc30(1,1)==0
                ;
            else
                % Plot the intercept point on boundary
                Nd1=NODE{fc33(1,2)};
                Nd2=NODE{fc33(2,2)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc33(3,2)*L*v;
                plot(C_i(1),C_i(2),'Marker', 'd','Markersize',3, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
                pause(1)
                hold on;
            end
        else
            c_d=CELL{fc30(1,2)};
            C_d=c_d{5};
            plot(C_d(1),C_d(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'g');
            pause(1)
            hold on;
        end
        plot(fc31(1,2),fc31(2,2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
        pause(1)
        hold on;
        %% Plot the upwind stencil point
        if fc30(1,3)==0
            % Plot the intercept point on boundary
            Nd1=NODE{fc33(1,3)};
            Nd2=NODE{fc33(2,3)};
            L=dis(Nd1{3},Nd2{3});
            v=Nd2{3}-Nd1{3};
            v=v/norm(v);
            C_i=Nd1{3}+fc33(3,3)*L*v;
            plot(C_i(1),C_i(2),'Marker', 'x','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
            pause(1)
            hold on;
        else
            c_u=CELL{fc30(1,3)};
            C_u=c_u{5};
            plot(C_u(1),C_u(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'm');
            pause(1)
            hold on;
        end
        plot(fc31(1,3),fc31(2,3),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm','MarkerEdgeColor', 'm');
        pause(1)
        hold on;
        %% Plot the further upwind stencil point
        if fc30(1,4)==0
            % Plot the intercept point on boundary
            if fc30(1,3)==0
                ;
            else
                Nd1=NODE{fc33(1,4)};
                Nd2=NODE{fc33(2,4)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc33(3,4)*L*v;
                plot(C_i(1),C_i(2),'Marker', 'x','Markersize',4, 'MarkerFaceColor', 'y','MarkerEdgeColor', 'y');
                pause(1)
                hold on;
            end
        else
            c_uu=CELL{fc30(1,4)};
            C_uu=c_uu{5};
            plot(C_uu(1),C_uu(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'r');
            pause(1)
            hold on;
        end
        plot(fc31(1,4),fc31(2,4),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        pause(1)
        hold on;
    end
end
