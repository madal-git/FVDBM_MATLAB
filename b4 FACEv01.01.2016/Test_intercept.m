N=360;
dt=2*pi/N;
n=zeros(2,N);
C=zeros(2,N);
ang=zeros(1,N);
L=2;
for i=1:N
    n(1,i)=cos((i-1)*dt);
    n(2,i)=sin((i-1)*dt);
    C(:,i)=[0;0]+L*n(:,i);
end

for i=1:N
    if i==N
        ang(i)=angle(n(:,i),n(:,1));
    else
        ang(i)=angle(n(:,i),n(:,i+1));
    end
end
%% Test angle
bar(ang)

Ori=[1;-2];
%% Test horizontal segment to be intercept
C1=[-14;10];
C2=[7.65;-16];

intercept_found_counter=0;
intercept_coord=zeros(2,1);

for i=1:N
    [Found, C_i] = intercept(Ori,n(:,i),C1,C2);
    if Found
        intercept_found_counter=intercept_found_counter+1;
        intercept_coord(:,intercept_found_counter)=C_i;
    end
end

figure
plot([C1(1),C2(1)],[C1(2),C2(2)],'black', 'linewidth',1);
hold on
plot(C1(1,1),C1(2,1),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'blue','MarkerEdgeColor', 'blue');
hold on
plot(C2(1,1),C2(2,1),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'blue','MarkerEdgeColor', 'blue');
hold on
plot(Ori(1,1),Ori(2,1),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'black','MarkerEdgeColor', 'black');
axis equal tight;
hold on
for s=1:intercept_found_counter;
    plot(intercept_coord(1,s),intercept_coord(2,s),'Marker', 'o','Markersize',2, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
    hold on;
end

