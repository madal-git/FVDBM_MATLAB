if FM==0
    T=0.5;
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
% l=1371;
for i=l:l
    FC=FACE{i};
    C_b=FC{7};
%     if fc23==1
    if FC{2}==0 % Check inner boundary
%             if FC{2}==1 % Check outter boundary
        %     if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose further downwind stencil point is out of boundary
        %     if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose further upwind stencil point is out of boundary
        %     if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose both further downwind  and further upwind stencil points are out of boundary
        %     if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose all stencil points are within the boundary
        %% Plot the stencil
        p1=C_b-FC{4}'*d;
        p2=C_b+FC{4}'*d;
        plot([p1(1),p2(1)],[p1(2),p2(2)],'blue', 'linewidth',1);
        pause(T)
        hold on
    end
end