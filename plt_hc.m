xlin=linspace(X1,X2,100);
ylin=linspace(Y1,Y2,100);
[Xx,Yy]=meshgrid(xlin,ylin);
if FTH==1
    figure;
    % use barycenter value
%     Z=griddata(XXX,YYY,T_plt,Xx,Yy);
%     contourf(Xx,Yy,Z,100);
    % use nodal value
    Z=griddata(X,Y,T_nd,Xx,Yy);
    contourf(Xx,Yy,Z,100);
    axis equal tight;
%     hold on
%     for r=1:M
%         P=CELL{r};
%         plot(P{22},P{23},'black',P{24},P{25},'black',P{26},P{27},'black','linewidth',0.5);
%         hold on
%     end
end
figure;
surf(Z)
figure;
plot(TR)