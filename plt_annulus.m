figure(4);
xlin=linspace(X1,X2,100);
ylin=linspace(Y1,Y2,100);
[Xx,Yy]=meshgrid(xlin,ylin);
Z=griddata(XXX,YYY,Rho_plt,Xx,Yy);
surf(Xx,Yy,Z);
axis equal tight;
hold on;
plot3(XXX,YYY,Rho_plt,'.','MarkerSize',15);
figure(5);
contourf(Xx,Yy,Z,100);
axis equal tight;
figure(6);
Z=griddata(XXX,YYY,U_M,Xx,Yy);
contourf(Xx,Yy,Z,100);
axis equal tight;
for r=1:M
    P=CELL{r};
    Centroid=P{5};
    XX(r)=Centroid(1,1);
    YY(r)=Centroid(2,1);
end
figure(10);
quiver(XX,YY,U_plt(1,1:M),U_plt(2,1:M),10);
axis equal tight  
figure(11);
plot(TT,log10(R));
figure(12);
plot(TT,RHO);
figure(13);
plot(TT,UR(1,:),TT,UR(2,:));
hold off


U_ref=Radius_outer*w_outer;
CX=Inf;
CY=(Y2+Y1)/2;
cut=[CX;CY];
[C_cut,U_cut]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,U,U_nd,CELL,M,NODE,N);
N_cut=length(U_cut);
% Inner boundary marker
Coord_inner_cylinder=1+2/6;
for i=1:N_cut
    if single(1000*epsilon+C_cut(1,i))==single(1000*epsilon+Coord_inner_cylinder)
        Marker_inner_cylinder=i;
        break;
    end
end
if Marker_inner_cylinder==N_cut
    error('Inner Boundary not found!');
end
% Outer boundary marker
Coord_outer_cylinder=1+Radius_outer;
for i=1:N_cut-1
    if C_cut(1,i)<Coord_outer_cylinder && C_cut(1,i+1)>Coord_outer_cylinder
        Marker_outer_cylinder=i;
        break;
    end
end
if Marker_outer_cylinder==N_cut-1
    error('Outer Boundary not found!');
end
% find the boundary velocity
bc_cell=CELL{in_which_triangle([Coord_outer_cylinder;1],CELL,M)};
bc_nd1=0;
bc_nd2=0;
bc_nd_counter=0;
for i=1:3
    ND=NODE{bc_cell{6+i}};
    if ND{2}==5
        bc_nd_counter=bc_nd_counter+1;
    end
    if bc_nd_counter==1
        bc_nd1=ND{1};
    elseif bc_nd_counter==2
        bc_nd2=ND{1};
    elseif bc_nd_counter==3
        error('logic error!');
    end
end
% velocity of slice
ND1=NODE{bc_nd1};
ND2=NODE{bc_nd2};
[Bool_found, C_i] = intercept([1;1],[1;0],ND1{3},ND2{3});
% coordinates of the slice
Coord_sim=[C_cut(1,Marker_inner_cylinder:Marker_outer_cylinder),C_i(1)];
figure
plot(C_i(1),C_i(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
hold on
Coord1=ND1{3};
plot(Coord1(1),Coord1(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'b','MarkerEdgeColor', 'b');
hold on
Coord2=ND2{3};
plot(Coord2(1),Coord2(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
hold on
axis equal tight
Dis_sum=dis(Coord1,Coord2);
Dis_seg1=dis(Coord1,C_i);
Dis_seg2=dis(Coord2,C_i);
if Dis_seg1+Dis_seg2~=Dis_sum
    error('Logic Error?');
end
c1=Dis_seg1/Dis_sum;
c2=Dis_seg2/Dis_sum;
U_bc_nd=c2*U_nd(:,ND1{1})+c1*U_nd(:,ND2{1});
U_sim=[U_cut(1,Marker_inner_cylinder:Marker_outer_cylinder),U_bc_nd(1,1)];
V_sim=[U_cut(2,Marker_inner_cylinder:Marker_outer_cylinder),U_bc_nd(2,1)];
U_ana=zeros(1,length(Coord_sim));
V_ana=zeros(1,length(Coord_sim));
for i=1:length(Coord_sim)
    V_ana(i)=((Coord_sim(i)-1)*(w_outer*(Coord_outer_cylinder-1)^2-w_inner*(Coord_inner_cylinder-1)^2)+(Coord_outer_cylinder-1)^2*(Coord_inner_cylinder-1)^2*(w_inner-w_outer)/(Coord_sim(i)-1))/((Coord_outer_cylinder-1)^2-(Coord_inner_cylinder-1)^2);
end

figure
plot((Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),V_sim,(Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),V_ana);
figure
plot((Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),V_sim-V_ana,(Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),U_sim-U_ana);
title('The velocity error compared to analytical solution')
xlabel('r/(R_o-R_i)')
ylabel('Velocity Error')
legend('x-velocity','y-velocity')
figure
plot((Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),(V_sim-V_ana)/norm(V_ana),(Coord_sim-Coord_sim(1))/(Coord_sim(end)-Coord_sim(1)),(U_sim-U_ana)/norm(V_ana));
title('Normalized velocity error compared to analytical solution')
xlabel('r/(R_o-R_i)')
ylabel('Normalized Velocity Error (u(v)/norm(v_ana))')
legend('x-velocity','y-velocity')

