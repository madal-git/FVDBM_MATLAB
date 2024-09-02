xlin=linspace(X1,X2,1000);
ylin=linspace(Y1,Y2,1000);
[Xx,Yy]=meshgrid(xlin,ylin);

x_2d=griddata(XXX,YYY,XXX,Xx,Yy);
y_2d=griddata(XXX,YYY,YYY,Xx,Yy);
u_2d=griddata(XXX,YYY,U_plt(1,:),Xx,Yy);
v_2d=griddata(XXX,YYY,U_plt(2,:),Xx,Yy);

figure
% quiver(x_2d,y_2d,u_2d,v_2d,10)
% axis equal tight  

% startx = 0:0.1:20;
% starty = ones(size(startx));
starty = 0:0.02:2;
startx = 0:0.02:2;

streamline(x_2d,y_2d,u_2d,v_2d,startx,starty)
axis equal tight  

hold on
zoomcenter(3,2.7,8)