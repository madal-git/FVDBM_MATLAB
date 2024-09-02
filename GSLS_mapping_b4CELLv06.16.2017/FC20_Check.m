e=1000;
a=0;
b=0;
for l=1:O
    FC=FACE{l};
    Map=FC{20};
    Map_down=Map(:,1);
    Map_up=Map(:,2);
    if Map_up(1,1)==1
        a=a+1;
    end
    if Map_down(1,1)==1
        b=b+1;
    end
end

figure(1);
for r=1:M;
    P=CELL{r};
    plot(P{22},P{23},P{24},P{25},P{26},P{27});
    hold on
    Cent=P{5};
    plot(Cent(1,1),Cent(2,1), 'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
    hold on
end;
axis equal tight;


for l=1:O
    FC=FACE{l};
    if FC{2}==0
        neigh_down=FC{13};
        Cell_down=CELL{neigh_down};
        C_down_fixed=Cell_down{5};
        neigh_up=FC{12};
        Cell_up=CELL{neigh_up};
        C_up_fixed=Cell_up{5};
        C=FC{18};
        C_down=C(:,1);
        C_up=C(:,2);
        plot(C_down(1,1),C_down(2,1), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
        hold on
        plot(C_up(1,1),C_up(2,1), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'blue','MarkerEdgeColor', 'blue');
        hold on
        Map=FC{20};
        Map_down=Map(:,1);
        Map_up=Map(:,2);
        %% Down
        if Map_down(1,1)==0
            if (Map_down(2,1)==Map_down(3,1)) && (Map_down(2,1)==Map_down(4,1))
                if single(e+norm(C_down-C_down_fixed))~=single(e)
                    error('The stencil point is supposed to be located at the centroid!');
                end
            else
                Cell1=CELL{Map_down(2,1)};
                Cell2=CELL{Map_down(3,1)};
                Cell3=CELL{Map_down(4,1)};
                C1=Cell1{5};
                C2=Cell2{5};
                C3=Cell3{5};
                plot([C1(1,1);C2(1,1)],[C1(2,1);C2(2,1)],'red', 'linewidth',1);
                hold on;
                plot([C2(1,1);C3(1,1)],[C2(2,1);C3(2,1)],'red', 'linewidth',1);
                hold on;
                plot([C3(1,1);C1(1,1)],[C3(2,1);C1(2,1)],'red', 'linewidth',1);
                hold on;
            end
        elseif Map_down(1,1)==1
            ND=NODE{Map_down(2,1)};
            Cell2=CELL{Map_down(3,1)};
            Cell3=CELL{Map_down(4,1)};
            C1=ND{3};
            C2=Cell2{5};
            C3=Cell3{5};
            plot([C1(1,1);C2(1,1)],[C1(2,1);C2(2,1)],'red', 'linewidth',1);
            hold on;
            plot([C2(1,1);C3(1,1)],[C2(2,1);C3(2,1)],'red', 'linewidth',1);
            hold on;
            plot([C3(1,1);C1(1,1)],[C3(2,1);C1(2,1)],'red', 'linewidth',1);
            hold on;
        elseif Map_down(1,1)==2
            Zone=FC{19};
            Zone_down=Zone(1,1);
            if Zone_down==1
                if Cell_down{10}==0 || Cell_down{11}==0
                    error('?');
                end
                if ~in_triangle(C_down,C_down_fixed,Cell_down{13},Cell_down{14})
                    error('?');
                end
            elseif Zone_down==2
                if Cell_down{11}==0 || Cell_down{12}==0
                    error('?');
                end
                if ~in_triangle(C_down,C_down_fixed,Cell_down{14},Cell_down{15})
                    error('?');
                end
            else
                if Cell_down{12}==0 || Cell_down{10}==0
                    error('?');
                end
                if ~in_triangle(C_down,C_down_fixed,Cell_down{15},Cell_down{13})
                    error('?');
                end
            end
        else
            error(':');
        end
        %% up
        if Map_up(1,1)==0
            if (Map_up(2,1)==Map_up(3,1)) && (Map_up(2,1)==Map_up(4,1))
                if single(e+norm(C_up-C_up_fixed))~=single(e)
                    error('The stencil point is supposed to be located at the centroid!');
                end
            else
                Cell1=CELL{Map_up(2,1)};
                Cell2=CELL{Map_up(3,1)};
                Cell3=CELL{Map_up(4,1)};
                C1=Cell1{5};
                C2=Cell2{5};
                C3=Cell3{5};
                plot([C1(1,1);C2(1,1)],[C1(2,1);C2(2,1)],'blue', 'linewidth',1);
                hold on;
                plot([C2(1,1);C3(1,1)],[C2(2,1);C3(2,1)],'blue', 'linewidth',1);
                hold on;
                plot([C3(1,1);C1(1,1)],[C3(2,1);C1(2,1)],'blue', 'linewidth',1);
                hold on;
            end
        elseif Map_up(1,1)==1
            ND=NODE{Map_up(2,1)};
            Cell2=CELL{Map_up(3,1)};
            Cell3=CELL{Map_up(4,1)};
            C1=ND{3};
            C2=Cell2{5};
            C3=Cell3{5};
            plot([C1(1,1);C2(1,1)],[C1(2,1);C2(2,1)],'blue', 'linewidth',1);
            hold on;
            plot([C2(1,1);C3(1,1)],[C2(2,1);C3(2,1)],'blue', 'linewidth',1);
            hold on;
            plot([C3(1,1);C1(1,1)],[C3(2,1);C1(2,1)],'blue', 'linewidth',1);
            hold on;
        elseif Map_up(1,1)==2
            Zone=FC{19};
            Zone_up=Zone(1,2);
            if Zone_up==1
                if Cell_up{10}==0 || Cell_up{11}==0
                    error('?');
                end
                if ~in_triangle(C_up,C_up_fixed,Cell_up{13},Cell_up{14})
                    error('?');
                end
            elseif Zone_up==2
                if Cell_up{11}==0 || Cell_up{12}==0
                    error('?');
                end
                if ~in_triangle(C_up,C_up_fixed,Cell_up{14},Cell_up{15})
                    error('?');
                end
            else
                if Cell_up{12}==0 || Cell_up{10}==0
                    error('?');
                end
                if ~in_triangle(C_up,C_up_fixed,Cell_up{15},Cell_up{13})
                    error('?');
                end
            end
        else
            error(':');
        end
    end
end
%