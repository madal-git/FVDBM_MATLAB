
%% Method 1
% temp=0;
% for i=1:N
%     temp=temp+U_nd(2,i)*T_nd(1,i);
% end
% temp=temp/N;
% Nu1=1+temp/(Tau_t/3*2*abs(T_bc(1)-T_bc(3))/(Y2-Y1))

%% Method 2, Calculate the Nu on the top and bottom surfaces
% Find the cells that are attached to the bottom surface and calculate the Nu on the bottom surface
e=10;
k_bottom=0;
BC_cell_bottom=0;
for i=1:M
    CL=CELL{i};
    for j=1:3
        ND=NODE{CL{6+j}};
        Coor=ND{3};
        if single(dis(Coor,[Coor(1,1);Y1])+e)==single(e)
            k_bottom=k_bottom+1;
            BC_cell_bottom(k_bottom)=i;
            break;
        end
    end
end

% figure;
% for r=1:k_bottom
%     P=CELL{BC_cell_bottom(r)};
%         plot(P{22},P{23},P{24},P{25},P{26},P{27});
%         hold on
% end
% axis equal tight

Nu_bottom=zeros(1,k_bottom);
T_s=T_bc(3);
T_inf=T_bc(1);
for i=1:k_bottom
    CL=CELL{BC_cell_bottom(i)};
    Coor=CL{5};
%     dis(Coor,[Coor(1,1);Y1])
    Nu_bottom(i)=((T_s-T(BC_cell_bottom(i)))/dis(Coor,[Coor(1,1);Y1]))/((T_s-T_inf)/(Y2-Y1));
end
Nu_bottom=sum(Nu_bottom)/k_bottom
%% Find the cells that are attached to the top surface and calculate the Nu on the top surface
k_top=0;
BC_cell_top=0;
for i=1:M
    CL=CELL{i};
    for j=1:3
        ND=NODE{CL{6+j}};
        Coor=ND{3};
        if single(dis(Coor,[Coor(1,1);Y2])+e)==single(e)
            k_top=k_top+1;
            BC_cell_top(k_top)=i;
            break;
        end
    end
end

% figure;
% for r=1:k_top
%     P=CELL{BC_cell_top(r)};
%         plot(P{22},P{23},P{24},P{25},P{26},P{27});
%         hold on
% end
% axis equal tight

Nu_top=zeros(1,k_top);
T_s=T_bc(3);
T_inf=T_bc(1);
for i=1:k_top
    CL=CELL{BC_cell_top(i)};
    Coor=CL{5};
%     dis(Coor,[Coor(1,1);Y2])
    Nu_top(i)=((T(BC_cell_top(i))-T_inf)/dis(Coor,[Coor(1,1);Y2]))/((T_s-T_inf)/(Y2-Y1));
end
Nu_top=sum(Nu_top)/k_top
Nu_avg=(Nu_top+Nu_bottom)/2
%% Method 3, calculate the Nu at each point, then calculate the average
% T_simu=zeros(N_H,N_L);
% if FM==0 %IRT mesh
%     if mod(N_L,2)==1
%         for j=1:N_L
%             C_cut_y=zeros(2,N_H);
%             C_cut_y(1,:)=X1+(X2-X1)/(N_L-1)*(j-1);
%             C_cut_y(2,:)=Y1:(Y2-Y1)/(N_H-1):Y2;
%             for i=1:N_H
%                 for l=1:N
%                     ND=NODE{l};
%                     if single(10+dis(C_cut_y(:,i),ND{3}))==single(10)
%                         T_simu(i,j)=T_nd(1,l);
%                         break;
%                     end
%                 end
%             end
%         end
%     else
%         error('Temporarily not available!');
%     end
% elseif FM==1 % General Mesh
%     for i=1:N_H
%         if i==1
%             T_simu(N_H,:)=[T_nd(1,N),T_nd(1,N_I+1:N_I+N_L-1)];
%         elseif i==N_H
%             T_simu(1,:)=fliplr(T_nd(1,N-(N_H-1)-(N_L-1):N-(N_H-1)));
%         else
%             cut=[Inf;Y2-dy*(i-1)];
%             [C_cut_x,T_simu(N_H-(i-1),:)]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,T,T_nd,CELL,M,NODE,N);
%         end
%     end
% else
%     error('Wrong flag for mesh type!');
% end
% 
% % Calculate the Nu at each point
% T_bottom=T_bc(3);
% T_top=T_bc(1);
% Nu_simu=zeros(N_H,N_L);
% dy=(Y2-Y1)/(N_H-1);
% for j=1:N_L
%     for i=1:N_H
%         if i==1 % Bottom surface
%             Nu_simu(i,j)=((T_bottom-T_simu(i+1,j))/dy)/((T_bottom-T_top)/(Y2-Y1));
%         elseif i==N_H % Top surface
%             Nu_simu(i,j)=((T_simu(i-1,j)-T_top)/dy)/((T_bottom-T_top)/(Y2-Y1));
%         else
%             Nu_simu(i,j)=((T_simu(i-1,j)-T_simu(i+1,j))/dy/2)/((T_bottom-T_top)/(Y2-Y1));
%         end
%     end
% end
% surf(Nu_simu)
% Nu_avg=sum(sum(Nu_simu))/N_L/N_H