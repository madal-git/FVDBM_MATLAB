function [Nu, Angle]=plt_nc(CELL,NODE,M,N,X1,X2,Y1,Y2,U,U_nd,T,T_nd,T_bc,Body_force,Tau,Tau_t)

% This function numerically calculates the Nu number distribution on the
% circular cylinder surface from 0 degree to 180 degrees, where 0 degree is
% the bottom point and 180 is the top point after moving counterclock-wise

Diameter=0.6;

%% Calculate and plot Nu distribution on the cylinder surface
Nu_nd=0; % Container for the node number that will be used for Nu number calculation
k=0;
for r=1:M
    P=CELL{r};
    if P{37}==3
        for i=1:3
            ND=NODE{P{i+6}};
            if ND{2}==0
                k=k+1;
                Nu_nd(k)=ND{1};
            end
        end
    end
end

R_nd=zeros(1,k);
Nu_nd=sort(Nu_nd);
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
    R_nd(i)=dis(nd_coor,[(X2-X1)/2;(Y2-Y1)/2]);
end

R_nd=sum(R_nd)/k;
Nd_start=0;
l=0;
e=2;
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
%     (nd_coor(2)-(Y2-Y1)/2)/R_nd
%     asin((nd_coor(2)-(Y2-Y1)/2)/R_nd)/pi*180
    if abs(asin((nd_coor(2)-(Y2-Y1)/2)/R_nd)/pi*180+90)<e
        l=l+1;
        Nd_start(l)=i;
    end
end
if length(Nd_start)~=2
    error('adjust e');
end

Nu_nd=[Nu_nd(Nd_start(2):end),Nu_nd(1:Nd_start(1))];
if length(Nu_nd)~=k
    error('logic error!');
end

Angle=zeros(1,k);
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
    Angle(i)=asin((nd_coor(2)-(Y2-Y1)/2)/R_nd)/pi*180+90;
end

Nu_node_0to180=Nu_nd(1:floor(k/2));
Angle=Angle(1:floor(k/2));
p=length(Nu_node_0to180);

% figure(1);
% for r=1:M
%     P=CELL{r};
%     if dis(P{5},[(X2-X1)/2;(Y2-Y1)/2])<Diameter/2+0.3
%         plot(P{22},P{23},P{24},P{25},P{26},P{27});
%         hold on
%     end
% end
% axis equal tight
% hold on
% for i=1:p
%     ND=NODE{Nu_node_0to180(i)};
%     nd_coor=ND{3};
%     plot(nd_coor(1),nd_coor(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
% %     pause(0.1)
%     hold on
% end

Nu=zeros(1,p);
T_s=T_bc(5);
T_inf=T_bc(1);
for i=1:p
    Nu(i)=((T_s-T_nd(Nu_node_0to180(i)))/(R_nd-Diameter/2))/((T_s-T_inf)/Diameter);
end
Ra=abs(Body_force(2))*(1/T_bc(1))*(T_bc(5)-T_bc(1))*Diameter^3/(Tau/3)/(Tau_t/3);
Pr=(Tau/3)/(Tau_t/3);
figure(2)
plot(Angle,Nu/Ra^0.25)
title(['Nu number distribution on the cynlinder surface from 0 to 180 degrees at Ra=', num2str(Ra), ' and Pr=', num2str(Pr)]);
xlabel('\theta,^{\circ}');
ylabel('Nu/Ra^1^/^4');

%% Create a group of points on a horizontal line starting from the node on the cylinder surface with 90 degrees and plot variables along the line
nd_90_start=0;
for i=1:N
    ND=NODE{i};
    if ND{2}<0
        nd_coor=ND{3};
        ang1=asin((nd_coor(2)-(Y2+Y1)/2)/(Diameter/2))/pi*180;
        if single(e+ang1)==single(e)
            nd_90_start=i;
            break;
        end
    end
end
if nd_90_start==0
    error('logoc error!');
end
c_b=nd_coor;
n_90=[1;0];
X_90_max=8;
S_coor=[0;0];
S_cl_id=0;
counter=0;
for i=1:M
    P=CELL{i};
    cl_coor=P{5};
    if abs(cl_coor(2)-(Y2+Y1)/2)<0.5
        c_j=norm_joint(c_b,n_90,cl_coor);
        if in_triangle(c_j,P{13},P{14},P{15})
            if c_j(1)<=X_90_max*Diameter/Ra^0.25+c_b(1) && c_j(1)>=c_b(1)
                counter=counter+1;
                S_coor(:,counter)=c_j;
                S_cl_id(counter)=i;
            end
        end
    end
end
S_coor_sort=zeros(2,counter);
S_cl_id_sort=zeros(1,counter);
x_sort=sort(S_coor(1,:));
for i=1:counter
    x=x_sort(i);
    for j=1:counter
        if single(S_coor(1,j))==single(x)
            S_coor_sort(:,i)=S_coor(:,j);
            S_cl_id_sort(i)=S_cl_id(j);
            break;
        end
    end
end

% figure(1)
% for i=1:counter
%     P=CELL{S_cl_id_sort(i)};
%     coord=P{5};
%     plot(coord(1),coord(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
%     pause(1)
%     hold on
%     plot(S_coor_sort(1,i),S_coor_sort(2,i), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
%     pause(1)
%     hold on
% end

X=[0,(S_coor_sort(1,:)-c_b(1))/Diameter*Ra^0.25];
L=length(X);
v=zeros(1,L);
t=zeros(1,L);
for i=1:L
    if i==1
        v(i)=0;
        t(i)=T_bc(5);
    else
        v(i)=U(2,S_cl_id_sort(i-1));
        t(i)=T(1,S_cl_id_sort(i-1));
    end
end

% Plot tangential velocity
figure(3)
plot(X,v*Diameter/Ra^0.5/(Tau_t/3))
title(['Angular velocity on 90-degree line from the surface at Ra=', num2str(Ra), ' and Pr=', num2str(Pr)]);
xlabel('X^*');
ylabel('V^*');

% Plot temperature
figure(4)
plot(X,(t-T_bc(1))/(T_bc(5)-T_bc(1)))
title(['Dimensionless temperature on the 90-degree line away from the surface at Ra=', num2str(Ra), ' and Pr=', num2str(Pr)]);
xlabel('X^*');
ylabel('T^*');

%% Create a group of points on a vertical line starting from the node on the cylinder surface with 180 degrees and plot variables along the line
nd_180_start=0;
for i=1:N
    ND=NODE{i};
    if ND{2}<0
        nd_coor=ND{3};
        ang1=abs(asin((nd_coor(2)-(Y2+Y1)/2)/(Diameter/2)))/pi*180;
        if single(e+ang1)==single(e+90)
            nd_180_start=i;
            break;
        end
    end
end
if nd_180_start==0
    error('logoc error!');
end
c_b=nd_coor;
n_180=[0;1];
Y_180_max=8;
S_coor_180=[0;0];
S_cl_id_180=0;
counter=0;
for i=1:M
    P=CELL{i};
    cl_coor=P{5};
    if abs(cl_coor(1)-(X2+X1)/2)<0.5
        c_j=norm_joint(c_b,n_180,cl_coor);
        if in_triangle(c_j,P{13},P{14},P{15})
            if c_j(2)<=Y_180_max*Diameter/Ra^0.25+c_b(2) && c_j(2)>=c_b(2)
                counter=counter+1;
                S_coor_180(:,counter)=c_j;
                S_cl_id_180(counter)=i;
            end
        end
    end
end
S_coor_180_sort=zeros(2,counter);
S_cl_id_180_sort=zeros(1,counter);
y_sort=sort(S_coor_180(2,:));
for i=1:counter
    y=y_sort(i);
    for j=1:counter
        if single(S_coor_180(2,j))==single(y)
            S_coor_180_sort(:,i)=S_coor_180(:,j);
            S_cl_id_180_sort(i)=S_cl_id_180(j);
            break;
        end
    end
end

% figure(1)
% for i=1:counter
%     P=CELL{S_cl_id_180_sort(i)};
%     coord=P{5};
%     plot(coord(1),coord(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
%     pause(1)
%     hold on
%     plot(S_coor_180_sort(1,i),S_coor_180_sort(2,i), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'g','MarkerEdgeColor', 'g');
%     pause(1)
%     hold on
% end

Y=[0,(S_coor_180_sort(2,:)-c_b(2))/Diameter*Ra^0.25];
L=length(Y);
v_180=zeros(1,L);
t_180=zeros(1,L);
for i=1:L
    if i==1
        v_180(i)=0;
        t_180(i)=T_bc(5);
    else
        v_180(i)=U(2,S_cl_id_180_sort(i-1));
        t_180(i)=T(1,S_cl_id_180_sort(i-1));
    end
end
figure(5)
plot(Y,v_180*Diameter/Ra^0.25/(Tau_t/3))
title(['Radial velocity on 180-degree line from the surface at Ra=', num2str(Ra), ' and Pr=', num2str(Pr)]);
xlabel('Y^*');
ylabel('U^*');

% Plot temperature
figure(4)
hold on
plot(Y,(t_180-T_bc(1))/(T_bc(5)-T_bc(1)))