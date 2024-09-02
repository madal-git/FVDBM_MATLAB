function paraview_streamline(U_0,N_L,N_H,N_I,X1,X2,Y1,Y2,U,U_nd,CELL,M,NODE,N)

u=zeros(N_H,N_L);
v=u;
X=zeros(N_H,N_L);
Y=X;
dx=(X2-X1)/(N_L-1);
dy=(Y2-Y1)/(N_H-1);
for i=1:N_H
    if i==1
        u(1,:)=[U_nd(1,N),U_nd(1,N_I+1:N_I+N_L-1)]/U_0;
        v(1,:)=[U_nd(2,N),U_nd(2,N_I+1:N_I+N_L-1)]/U_0;
        X(1,:)=(X1:dx:X2)/(X2-X1);
        Y(1,:)=ones(1,N_L)*Y2/(Y2-Y1);
    elseif i==N_H
        u(N_H,:)=fliplr(U_nd(1,N-(N_H-1)-(N_L-1):N-(N_H-1)))/U_0;
        v(N_H,:)=fliplr(U_nd(2,N-(N_H-1)-(N_L-1):N-(N_H-1)))/U_0;
        X(N_H,:)=(X1:dx:X2)/(X2-X1);
        Y(N_H,:)=ones(1,N_L)*Y1/(Y2-Y1);
    else
        cut=[Inf;Y2-dy*(i-1)];
        [C_cut_x,U_cut]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,U,U_nd,CELL,M,NODE,N);
        u(i,:)=U_cut(1,:)/U_0;
        v(i,:)=U_cut(2,:)/U_0;
        X(i,:)=C_cut_x(1,:)/(X2-X1);
        Y(i,:)=C_cut_x(2,:)/(Y2-Y1);
    end
end
quiver(X,Y,u,v,10);
axis equal tight

save('Cavity_streamline.mat','u','v','X','Y');

B=cat(3,u,v);
B=permute(B,[3 1 2]);
fid = fopen('paraview_streamline.raw','w');
fwrite(fid,B,'float'); 
fclose(fid);