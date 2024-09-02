xlin=linspace(X1,X2,100);
ylin=linspace(Y1,Y2,100);
[Xx,Yy]=meshgrid(xlin,ylin);
if FTH==1
    figure(1)
    % use barycenter value
%     Z=griddata(XXX,YYY,T_plt,Xx,Yy);
%     contourf(Xx,Yy,Z,100);
    % use nodal value
    Z=griddata(X,Y,T_nd,Xx,Yy);
    contourf(Xx,Yy,Z,100);
    axis equal tight;
end
figure
surf(Z)