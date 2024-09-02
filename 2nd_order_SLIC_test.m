a=0.1;
b=0.05;
for r=1:M
    Cell=CELL{r};
    coor=Cell{5};
    f_old(1:qh,r)=[0.5;0.1;0.2;0.3;0.4;0.01;0.02;0.03;0.04]+(a*coor(1,1)+b*coor(2,1));
    f_im_ana(:,r)=f_old(1:qh,r);
end

f_im=in_cell_mapping_SLIC_apply(CELL,V,f_old,FPDC);


XXX=zeros(1,M);
YYY=zeros(1,M);
%%%% coordinates
for r=1:M
    P=CELL{r};
    Centroid=P{5};
    XXX(r)=Centroid(1,1);
    YYY(r)=Centroid(2,1);
end
for i=1:qh
    figure;
    xlin=linspace(X1,X2,100);
    ylin=linspace(Y1,Y2,100);
    [Xx,Yy]=meshgrid(xlin,ylin);
    Z=griddata(XXX,YYY,f_im(i,:),Xx,Yy);
    contourf(Xx,Yy,Z,100);
    axis equal tight;
    figure;
    Z=griddata(XXX,YYY,f_im(i,:)-f_im_ana(i,:),Xx,Yy);
    contourf(Xx,Yy,Z,100);
    axis equal tight;
end
