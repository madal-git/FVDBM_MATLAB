function plt_cl(X1,X2,Y1,Y2,CELL,M,NODE)
figure;

plot([X1,X2],[Y2,Y2],'black', 'linewidth',1);
hold on
plot([X2,X2],[Y2,Y1],'black', 'linewidth',1);
hold on
plot([X2,X1],[Y1,Y1],'black', 'linewidth',1);
hold on
plot([X1,X1],[Y1,Y2],'black', 'linewidth',1);
hold on
axis equal tight

for r=1:M
    CL=CELL{r};
    if CL{37}==2
        cell_neigh=CL{44};
        on_bc_counter=0;
        for i=1:3
            if cell_neigh(i)==0
                on_bc_counter=on_bc_counter+1;
            end
        end
        if on_bc_counter==0
            ;
        elseif on_bc_counter==1
            nd=CL{45};
            point1=CL{5};
            nd1=NODE{nd(1)};
            nd2=NODE{nd(2)};
            point2=nd1{3};
            point3=nd2{3};
            plot([point1(1),point2(1)],[point1(2),point2(2)],'red', 'linewidth',1);
            hold on
            pause(0.2)
            plot([point2(1),point3(1)],[point2(2),point3(2)],'blue', 'linewidth',1);
            hold on
            pause(0.2)
            plot([point3(1),point1(1)],[point3(2),point1(2)],'green', 'linewidth',1);
            hold on
            pause(0.2)
        else
            error('Logic error');
        end
    end
end

% for r=1:M
%     CL=CELL{r};
%     if CL{37}==2
%         cell_neigh=CL{44};
%         on_bc_counter=0;
%         for i=1:3
%             if cell_neigh(i)==0
%                 on_bc_counter=on_bc_counter+1;
%             end
%         end
%         if on_bc_counter==0
%             ;
%         elseif on_bc_counter==1
%             nd=CL{45};
%             point1=CL{5};
%             nd1=NODE{nd(1)};
%             nd2=NODE{nd(2)};
%             point2=nd1{3};
%             point3=nd2{3};
%             plot([point1(1),point2(1)],[point1(2),point2(2)],'red', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point2(1),point3(1)],[point2(2),point3(2)],'blue', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point3(1),point1(1)],[point3(2),point1(2)],'green', 'linewidth',1);
%             hold on
%             pause(0.2)
%         else
%             error('Logic error');
%         end
%     end
% end


% for r=1:M
%     CL=CELL{r};
%     if CL{37}==2
%         cell_neigh=CL{58};
%         on_bc_counter=0;
%         for i=1:3
%             if cell_neigh(i)==0
%                 on_bc_counter=on_bc_counter+1;
%             end
%         end
%         if on_bc_counter==0
%             ;
%         elseif on_bc_counter==1
%             nd=CL{59};
%             point1=CL{5};
%             nd1=NODE{nd(1)};
%             nd2=NODE{nd(2)};
%             point2=nd1{3};
%             point3=nd2{3};
%             plot([point1(1),point2(1)],[point1(2),point2(2)],'red', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point2(1),point3(1)],[point2(2),point3(2)],'blue', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point3(1),point1(1)],[point3(2),point1(2)],'green', 'linewidth',1);
%             hold on
%             pause(0.2)
%         else
%             error('Logic error');
%         end
%     end
% end

% for r=1:M
%     CL=CELL{r};
%     if CL{37}==2
%         cell_neigh=CL{65};
%         on_bc_counter=0;
%         for i=1:3
%             if cell_neigh(i)==0
%                 on_bc_counter=on_bc_counter+1;
%             end
%         end
%         if on_bc_counter==0
%             ;
%         elseif on_bc_counter==1
%             nd=CL{66};
%             point1=CL{5};
%             nd1=NODE{nd(1)};
%             nd2=NODE{nd(2)};
%             point2=nd1{3};
%             point3=nd2{3};
%             plot([point1(1),point2(1)],[point1(2),point2(2)],'red', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point2(1),point3(1)],[point2(2),point3(2)],'blue', 'linewidth',1);
%             hold on
%             pause(0.2)
%             plot([point3(1),point1(1)],[point3(2),point1(2)],'green', 'linewidth',1);
%             hold on
%             pause(0.2)
%         else
%             error('Logic error');
%         end
%     end
% end