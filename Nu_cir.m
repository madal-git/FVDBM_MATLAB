Diameter=0.6;
figure;
for r=1:M
    P=CELL{r};
    if dis(P{5},[(X2-X1)/2;(Y2-Y1)/2])<Diameter+0.04
        plot(P{22},P{23},P{24},P{25},P{26},P{27});
        hold on
    end
end
axis equal tight

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
figure
hold on
R_nd=zeros(1,k);
Nu_nd=sort(Nu_nd);
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
    R_nd(i)=dis(nd_coor,[(X2-X1)/2;(Y2-Y1)/2]);
    plot(nd_coor(1),nd_coor(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
    hold on
    pause(0.1)
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
figure
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
    plot(nd_coor(1),nd_coor(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
    hold on
    pause(0.1)
end


Angle=zeros(1,k);
for i=1:k
    ND=NODE{Nu_nd(i)};
    nd_coor=ND{3};
    Angle(i)=asin((nd_coor(2)-(Y2-Y1)/2)/R_nd)/pi*180+90;
end

Nu_node_0to180=Nu_nd(1:floor(k/2));
Angle_0to180=Angle(1:floor(k/2));

figure
p=length(Nu_node_0to180);
for i=1:p
    ND=NODE{Nu_node_0to180(i)};
    nd_coor=ND{3};
    plot(nd_coor(1),nd_coor(2), 'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
    hold on
    pause(0.1)
end


Nu=zeros(1,p);
T_s=T_bc(5);
T_inf=T_bc(1);
for i=1:p
    Nu(i)=((T_s-T_nd(Nu_node_0to180(i)))/(R_nd-Diameter/2))/((T_s-T_inf)/Diameter);
end
Ra=abs(Body_force(2))*Diameter^3/(Tau/3)/(Tau_t/3);
plot(Angle_0to180,Nu/Ra^0.25)