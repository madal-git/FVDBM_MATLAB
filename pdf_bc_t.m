function g_bc=pdf_bc_t(N_L,N_H,N_I,X1,X2,Y1,Y2,tau_g,dt,CELL,ND,f,feq,f_node,g,geq,g_node,U_node,U,RHO_node,RHO,T,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,qt,wc,wh,wt,fd,fdcy,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80)
% g_bc=pdf_bc_t(N_L,N_H,N_I,tau_g,dt,NODE,ND,g,geq,g_node,U_node,U,RHO_node,RHO,V,V1,V2,Rho_r,Rho_in,Rho_out,qt,wc,wh,fd,fdcy,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77) calculates
% the thermal pdf of the current node that is on boundary.
% N_I is the total number of interior nodes
% NODE is the entire NODE data structure
% ND is the NODE data structure of current node
% g is the thermal pdf matrix of triangle centroids of the entire domain
% geq is the thermal equilibrium pdf matrix of triangle centroids of the entire domain
% g_node is the pdf matrix of nodes of the entire domain
% U_node is the velocity matrix of nodes of the entire domain
% U is the velocity matrix of triangle centroid of the entire domain
% RHO is the density matrix of triangle centroid of the entire domain
% V is the lattice applied
% V1 is the first lattice structure
% V2 is the second lattice structure
% Rho_r is the reference flow density
% Rho_in is the flow density at inlet.
% Rho_out is the flow density at outlet.
% wc is weighting factor for calculating the nodal pdf at corner nodes by
% using extrapolation scheme
% For corner nodes have one side is periodic, wc is the weight of
% extrapolation in periodic direction
% wc*Extrapolation_periodic + (1-wc)*Extrapolation_nonperiodic
% For corner nodes without periodic on any side, wc is the weight of
% extrapolation on the upstream side.
% wc*Extrapolation_upstream + (1-wc)*Extrapolation_downstream
% wh is the hydrodynamic equilibrium weight coefficient for current lattice
% fd is the flag for which density is used. fd=0----local density; fd=1----reference density
% NPdc is the three digit scheme code for Periodic boundary condition.
% NInt is the four digit scheme code for Inlet boundary condition.
% NOut is the four digit scheme code for Outlet boundary condition.
% NStw is the three digit scheme code for Stationary Wall boundary condition.
% NMow is the three digit scheme code for Moving Wall boundary condition.
% NWdp is the three digit scheme code for Well Developed boundary condition.
% NIof is the three digit scheme code for Inflow/Outflow boundary condition.
% NC70 is the four digit scheme code for corner node with Periodic + Periodic boundary condition.
% NC71 is the four digit scheme code for corner node with Periodic + Velocity or Velocity + Periodic boundary condition.
% NC72 is the four digit scheme code for corner node with Periodic + Density or Density + Periodic boundary condition.
% NC73 is the four digit scheme code for corner node with Velocity + Velocity boundary condition.
% NC74 is the four digit scheme code for corner node with Density + Density boundary condition.
% NC75 is the four digit scheme code for corner node with Velocity + Density or Density + Velocity boundary condition.
% NC76 is the four digit scheme code for corner node with Well Developed + others or others + Well Developed boundary condition.
% NC77 is the four digit scheme code for corner node with Well Developed + Well Developed boundary condition.
% NC78 is the four digit scheme code for corner node with Periodic + Inflow/Outflow boundary condition.
% NC79 is the four digit scheme code for corner node with Density or Velocity + Inflow/Outflow boundary condition.
% NC80 is the four digit scheme code for corner node with Inflow/Outflow + Inflow/Outflow boundary condition.


R=8.314; % Thermal constant
q=length(g(:,1));
if q~=length(V(1,:)) || q~=length(wt)
    error('Check dimension of lattice, pdf and equilibrium coefficient!');
end
g_bc=zeros(q,1);

% e=ceil(10*norm(V)); %%%% Calculate the reference size


if ND{2}==0
    error('The current node is NOT on boundary!');
elseif abs(ND{2})==1 % Periodic
    if ND{2}<0
        error('Immersed boundary cannot be periodic!');
    end
    if NPdc==131 % local average
        g_bc=node_star(ND,g,1);
    else
        error('Input scheme code for nodal pdf on Periodic boundary is not valid or available!');
    end
elseif abs(ND{2})==20 % Velocity Inlet
    if ND{2}<0
        error('Immersed boundary cannot be velocity inlet!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    
    if NInt==2010 % local average Rho + local average f_neq
        %% Algo 1
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,0),T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    elseif NInt==2099 % Given U and Rho on boundary + local average f_neq
        g_bc_eq(:,1)=eqm_t(V,U_bc,Rho_bc,T_bc,q,wt,Rho_r,fd);
        g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
        g_bc=g_bc_eq+g_bc_neq;
    else
        error('Input scheme code for nodal pdf on Stationary Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==21 % Pressure Inlet
    if ND{2}<0
        error('Immersed boundary cannot be pressure inlet!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    T_bc=bc(4,1);
    if NInt==2100 % local average U + local average f_neq
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,node_star(ND,U,0),Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    else
        error('Input scheme code for nodal pdf on Density Inlet boundary is not valid or available!');
    end
elseif abs(ND{2})==30 % Velocity Outlet
    if ND{2}<0
        error('Immersed boundary cannot be velocity outlet!');
    end
elseif abs(ND{2})==31 % Pressure Outlet
    if ND{2}<0
        error('Immersed boundary cannot be pressure outlet!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    T_bc=bc(4,1);
    % f = f_eq + f_neq
    if NOut==3100 % local average U + local average f_neq
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,node_star(ND,U,0),Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    elseif NOut==3199 % v fixed to zero & extrapolated u + local average f_neq
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            U_ex=node_star(ND,U,0);
            g_bc_eq(:,1)=eqm_t(V,[U_ex(1);0],Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    end
elseif abs(ND{2})==4 % Stationary Wall
    bc=ND{21};
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    
    if NStw==410 % local average Rho + local average f_neq
%         if qh==qt
%             g_bc=R*T_bc*f_node(:,ND{1});
%         else
            g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,0),T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            % BB &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            if fdcy==0 % No decay
                g_bc=g_bc_eq+g_bc_neq;
            elseif fdcy==1 % with decay
                g_bc=g_bc_eq+g_bc_neq*exp(-dt/tau_g);
            else
                error('The flag for decaying of f_neq is invalid!');
            end
%         end
    elseif NStw==415 % f = f_eq + f_neq. Second-order mapping for the unknown part 
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_g=zeros(num_star_cell,q);
        F_geq=zeros(num_star_cell,q);
        for i=1:num_star_cell
            F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_g(i,:)=g(:,star_cell_index(i))';
            F_geq(i,:)=geq(:,star_cell_index(i))';
        end
        g_bc_eq(:,1)=eqm_t(V,U_bc,([1,ND{3}']*(CoeX*F_RHO))',T_bc,q,wt,Rho_r,fd);

        g_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_g))'-([1,ND{3}']*(CoeX*F_geq))';

        g_bc=g_bc_eq+g_bc_neq;
    else
        error('Input scheme code for nodal pdf on Stationary Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==5 % Moving Wall
    if ND{2}<0
        error('Immersed boundary cannot be moving wall!');
    end
    bc=ND{21};
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    
    if NMow==510 % local average Rho + local average f_neq
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,0),T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            if fdcy==0 % No decay
                g_bc=g_bc_eq+g_bc_neq;
            elseif fdcy==1 % with decay
                g_bc=g_bc_eq+g_bc_neq*exp(-dt/tau_g);
            else
                error('The flag for decaying of f_neq is invalid!');
            end
        end
    elseif NMow==515 % f = f_eq + f_neq. Second-order mapping for the unknown part 
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_g=zeros(num_star_cell,q);
        F_geq=zeros(num_star_cell,q);
        for i=1:num_star_cell
            F_RHO(i,:)=T(1,star_cell_index(i));
            F_g(i,:)=g(:,star_cell_index(i))';
            F_geq(i,:)=geq(:,star_cell_index(i))';
        end
        g_bc_eq(:,1)=eqm_t(V,U_bc,([1,ND{3}']*(CoeX*F_RHO))',T_bc,q,wt,Rho_r,fd);

        g_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_g))'-([1,ND{3}']*(CoeX*F_geq))';

        g_bc=g_bc_eq+g_bc_neq;
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==6 % Well Developed
    if ND{2}<0
        error('Immersed boundary cannot be Well Developed!');
    end
    % Zero gradient, f = f_upwind
    if NWdp==620
        g_bc=node_star(ND,g,0);
    elseif NWdp==640
        %% Determing the boundary outward normal;
        coord=ND{3};
        if single(coord(1,1))==single(X1) % Left boundary
            face_normal=[-1;0];
        elseif single(coord(1,1))==single(X2) % Right boundary
            face_normal=[1;0];
        elseif single(coord(2,1))==single(Y1) % Bottom boundary
            face_normal=[0;-1];
        elseif single(coord(2,1))==single(Y2) % Top boundary
            face_normal=[0;1];
        else
            error('Logic error');
        end
        %% Find the cell whose centroid is lined up with the face normal
        star_cell_count=ND{4};
        cell_star=ND{5};
        dot_product=zeros(1,star_cell_count);
        cell_marker=zeros(1,star_cell_count);
        for i=1:star_cell_count
            CL=CELL{cell_star(i)};
            dot_product(i)=(CL{5}-ND{3})'*face_normal;
        end
        min_dot_product=min(dot_product);
        for i=1:star_cell_count
            if single(10+dot_product(i))==single(10+min_dot_product)
                cell_marker(i)=1;
            end
        end
        if sum(cell_marker)~=1
            error('logic error!');
        end
        cell_found=cell_star*cell_marker';
        %% Calcuate the boundary condition
%         g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,1),T_bc,q,wt,Rho_r,fd);
%         g_bc_neq(:,1)=node_star(ND,g,1)-node_star(ND,geq,1);
%         g_bc=g_bc_eq+g_bc_neq;
        T_bc=T(1,cell_found);
        g_bc_eq(:,1)=eqm_t(V,[0;0],Rho_r,T_bc,q,wt,Rho_r,fd);
        g_bc_neq(:,1)=node_star(ND,g,1)-node_star(ND,geq,1);
        g_bc=g_bc_eq+g_bc_neq;
    else
        error('Input scheme code for nodal pdf on Well Developed boundary is not valid or available!');
    end
elseif abs(ND{2})==8 % Inflow/Outflow
    if ND{2}<0
        error('Immersed boundary cannot be Inflow/Outflow!');
    end
    if NIof==800
        g_bc=node_star(ND,g,0);
    else
        error('Input scheme code for nodal pdf on Inflow/Outflow boundary is not valid or available!');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Corner nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif abs(ND{2})==70 % Periodic + Periodic
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Periodic!');
    end
    if NC70==7000
        g_bc=node_star(ND,g,1);
    else
        error('Input scheme code for nodal pdf at two-sided periodic corner node is not valid or available!');
    end
elseif abs(ND{2})==71 % Periodic + Velocity or Velocity + Periodic
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Velocity or Velocity + Periodic!');
    end
    bc=ND{21};
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    
    if NC71==7110 % local average Rho + local average f_neq
%         if qh==qt
%             g_bc=R*T_bc*f_node(:,ND{1});
%         else
            g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,1),T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,1)-node_star(ND,geq,1);
            g_bc=g_bc_eq+g_bc_neq;
%         end
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==72 % Periodic + Density or Density + Periodic
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Density or Density + Periodic!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    T_bc=bc(4,1);
    
    if NC72==7210 % Local Rho + extrapolated U
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,node_star(ND,U,1),Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,1)-node_star(ND,geq,1);
            g_bc=g_bc_eq+g_bc_neq;
        end
    else
        error('Input scheme code for nodal pdf on periodic + density boundary is not valid or available!');
    end
elseif abs(ND{2})==73 % Velocity + Velocity
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    
    if NC73==7310 % local average Rho + local average f_neq
%         if qh==qt
%             g_bc=R*T_bc*f_node(:,ND{1});
%         else
            g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,0),T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            if fdcy==0 % No decay
                g_bc=g_bc_eq+g_bc_neq;
            elseif fdcy==1 % with decay
                g_bc=g_bc_eq+g_bc_neq*exp(-dt/tau_g);
            else
                error('The flag for decaying of f_neq is invalid!');
            end
%         end
    elseif NC73==7399 % both velocity and density are prescibed
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,U_bc,Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==74 % Density + Density
    if ND{2}<0
        error('Immersed boundary cannot be Density + Density!');
    end
elseif abs(ND{2})==75 % Velocity + Density or Density + Velocity
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Periodic!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    % f = f_eq + f_neq with weighted two-side extrapolation
    if NC75==7500 % local average f_neq
        if qh==qt
            g_bc=R*T_bc*f_node(:,ND{1});
        else
            g_bc_eq(:,1)=eqm_t(V,U_bc,Rho_bc,T_bc,q,wt,Rho_r,fd);
            g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
            g_bc=g_bc_eq+g_bc_neq;
        end
    else
        error('Input scheme code for nodal pdf on Velocity + Density or Density + Velocity corner node is not valid or available!');
    end
elseif abs(ND{2})==76 % Well Developed + others or others + Well Developed
    if ND{2}<0
        error('Immersed boundary cannot be Well Developed + others or others + Well Developed!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    % Zero gradient, f = f_upwind
    if NC76==7620
        g_bc=node_star(ND,g,0);
    elseif NC76==7621
        g_bc_eq(:,1)=eqm_t(V,U_bc,node_star(ND,RHO,0),node_star(ND,T,0),q,wt,Rho_r,fd);
        g_bc_neq(:,1)=node_star(ND,g,0)-node_star(ND,geq,0);
        % BB &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        if fdcy==0 % No decay
            g_bc=g_bc_eq+g_bc_neq;
        elseif fdcy==1 % with decay
            g_bc=g_bc_eq+g_bc_neq*exp(-dt/tau_g);
        else
            error('The flag for decaying of f_neq is invalid!');
        end
    else
        error('Input scheme code for nodal pdf on Well Developed + others corner node is not valid or available!');
    end
elseif abs(ND{2})==77 % Well Developed + Well Developed
    if NC77==7700
        g_bc=node_star(ND,g,0);
    else
        error('Input scheme code for nodal pdf on Well Developed + Well Developed corner node is not valid or available!');
    end
    if ND{2}<0
        error('Immersed boundary cannot be Well Developed + Well Developed!');
    end
elseif abs(ND{2})==78 % Periodic + Inflow/Outflow
    if NC78==7800
        g_bc=node_star(ND,g,1);
    else
        error('Check the flag of corner node!')
    end
elseif abs(ND{2})==79 % Velocity or Density + Inflow/Outflow
    if NC79==7900
        g_bc=node_star(ND,g,0);
    else
        error('Check the flag of corner node!')
    end
elseif abs(ND{2})==80 % Inflow/Outflow + Inflow/Outflow
    if NC80==8000
        g_bc=node_star(ND,g,0);
    else
        error('Check the flag of corner node!')
    end
else
    error('The boundary condition for current node is invalid or unavailable!');
end

% if B==1
%     for l=1:q
%         if V==V1
%             UB=ND{19};
%         elseif V==V2
%             UB=ND{20};
%         else
%             error('The lattice used is not the lattice used for determing upwind cells');
%         end
%         if f_bc_neq(l,1)==Inf
%             f_bc_neq(l,1)=f_bc_neq(UB(l),1); % replaced with correct number
%         end
%     end
%     f_bc=f_bc_eq+f_bc_neq;
% end