if FM==1 % General Mesh
    %% Numerical solution
    CX=Inf;
    CY=(Y2+Y1)/2;
    cut=[CX;CY];
    % Find the nodal velocity on the horizontal mid-slide
    [C_cut_x,U_cut]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,U,U_nd,CELL,M,NODE,N);
    % Selection the RHS section ranging from the inner boundary to ther outer boundary
    counter=0;
    mid_slice_x=0;
    U_mid_slice=[0;0];
    for i=1:length(C_cut_x)
        if C_cut_x(1,i)>Radius_inner+1 && C_cut_x(1,i)<Radius_outer+1
            counter=counter+1;
            mid_slice_x(counter)=C_cut_x(1,i);
            U_mid_slice(:,counter)=U_cut(:,i);
        end
    end
    % Fill the left and right end velocity on the boundary
    mid_slice_x=[Radius_inner+1,mid_slice_x,Radius_outer+1];
    U_mid_slice=[[0;Radius_inner*w_inner],U_mid_slice,[0;Radius_outer*w_outer]];
    % Analytical solution and plot
    r_sim=mid_slice_x-1;
    U_sim=U_mid_slice;
    u_ana=zeros(1,length(r_sim));
    v_ana=zeros(1,length(r_sim));
    for i=1:length(r_sim)
        v_ana(i)=(r_sim(i)*(w_outer*Radius_outer^2-w_inner*Radius_inner^2)+Radius_outer^2*Radius_inner^2*(w_inner-w_outer)/r_sim(i))/(Radius_outer^2-Radius_inner^2);
    end
end

figure;
plot((r_sim-Radius_inner)/(Radius_outer-Radius_inner),U_sim(2,:),(r_sim-Radius_inner)/(Radius_outer-Radius_inner),v_ana)
title('The Y-velocity comparison')
xlabel('r/(R_i_n_n_e_r-R_o_u_t_e_r)')
ylabel('Y-Velocity')
legend('FVDBM','Analytical')

figure
plot((r_sim-Radius_inner)/(Radius_outer-Radius_inner),U_sim(1,:),(r_sim-Radius_inner)/(Radius_outer-Radius_inner),u_ana)
title('The X-velocity comparison')
xlabel('r/(R_i_n_n_e_r-R_o_u_t_e_r)')
ylabel('X-Velocity')
legend('FVDBM','Analytical')

figure;
plot((r_sim-Radius_inner)/(Radius_outer-Radius_inner),(U_sim(1,:)-u_ana),(r_sim-Radius_inner)/(Radius_outer-Radius_inner),(U_sim(2,:)-v_ana))
title('The velocity error compared to analytical solution')
xlabel('r/(R_i_n_n_e_r-R_o_u_t_e_r)')
ylabel('Velocity Error')
legend('x-velocity','y-velocity')

figure;
plot((r_sim-Radius_inner)/(Radius_outer-Radius_inner),(U_sim(1,:)-u_ana)/norm(v_ana),(r_sim-Radius_inner)/(Radius_outer-Radius_inner),(U_sim(2,:)-v_ana)/norm(v_ana))
title('Normalized velocity error compared to analytical solution')
xlabel('r/(R_i_n_n_e_r-R_o_u_t_e_r)')
ylabel('Normalized Velocity Error (u(v)/norm(v_ana))')
legend('x-velocity','y-velocity')