function f_bc=pdf_bc_h(N_L,N_H,N_I,tau,dt,NODE,ND,CELL,f,feq,f_node,U_node,U,RHO_node,RHO,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,fd,fdcy,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80)
% f_bc=pdf_bc_h(ND,f,V,U_r,Rho_r,wh,NPdc,NInt,NOut,NStw,NMow) calculates
% the hydrodynamic pdf of the current node that is on boundary.
% N_I is the total number of interior nodes
% NODE is the entire NODE data structure
% ND is the NODE data structure of current node
% f is the pdf matrix of triangle centroids of the entire domain
% feq is the equilibrium pdf matrix of triangle centroids of the entire domain
% f_node is the pdf matrix of nodes of the entire domain
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

N=length(NODE);
q=length(f(:,1));
if q~=length(V(1,:)) || q~=length(wh)
    error('Check dimension of lattice, pdf and equilibrium coefficient!');
end
f_bc=zeros(q,1);

e=ceil(10*norm(V)); %%%% Calculate the reference size

B=0;
% B is Bounceback scheme flag
%%% B=0 ---- No bounceback scheme is applied
%%% B=1 ---- Bounceback scheme is applied

if ND{2}==0
    error('The current node is NOT on boundary!');
elseif abs(ND{2})==1 % Periodic
    if ND{2}<0
        error('Immersed boundary cannot be periodic!');
    end
    % f=f_eq
    if NPdc==100 % f=eqm(U,local average Rho)
        f_bc(:,1)=eqm_h(V,node_star(ND,U,1),node_star(ND,RHO,1),q,wh,Rho_r,fd);
    elseif NPdc==101 % f=eqm(U,linear extraplated Rho)
        NS=ND{16};
        NDL=NODE{NS(1)};
        NDR=NODE{NS(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,node_star(ND,U,1),0.5*(Rho_L+Rho_R),q,wh,Rho_r,fd);
        % f=f_eq + f_neq with flexible density
    elseif NPdc==110 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NPdc==111 % linear extrapolated Rho + local average f_neq
        NS=ND{16};
        NDL=NODE{NS(1)};
        NDR=NODE{NS(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),0.5*(Rho_L+Rho_R),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NPdc==112 % local average Rho + linear extrapolated f_neq
        NS=ND{16};
        NDL=NODE{NS(1)};
        NDR=NODE{NS(3)};
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=0.5*((f_node(:,NS(1))-f_eq_l(:,1))+(f_node(:,NS(3))-f_eq_r(:,1)));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NPdc==113 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{16};
        NDL=NODE{NS(1)};
        NDR=NODE{NS(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),0.5*(Rho_L+Rho_R),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=0.5*((f_node(:,NS(1))-f_eq_l(:,1))+(f_node(:,NS(3))-f_eq_r(:,1)));
        f_bc=f_bc_eq+f_bc_neq;
        % f=f_eq + f_neq with fixed density
    elseif NPdc==120 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NPdc==121 % linear extrapolated f_neq
        NS=ND{16};
        NDL=NODE{NS(1)};
        NDR=NODE{NS(3)};
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=0.5*((f_node(:,NS(1))-f_eq_l(:,1))+(f_node(:,NS(3))-f_eq_r(:,1)));
        f_bc=f_bc_eq+f_bc_neq;
        % f=f_total
    elseif NPdc==130 % copy of exterior neighbor
        NS=ND{16};
        f_bc(:,1)=f_node(:,NS(1));
    elseif NPdc==131 % local average
        f_bc=node_star(ND,f,1);
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
    
    % f = f_eq
    if NInt==2000 % fixed density
        f_bc(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
    elseif NInt==2001 % local average Rho
        f_bc(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
    elseif NInt==2002 % linear extrapolated Rho
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        % f = f_eq + f_neq
    elseif NInt==2010 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2011 % linear extrapolated Rho + local average f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2012 % local average Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2013 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with fixed density
    elseif NInt==2015 % Second-order
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,([1,ND{3}']*(CoeX*F_RHO))',T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2020 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2021 % linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + bounce-back f_neq
    elseif NInt==2030 % local average Rho + local average f_neq
        B=1; % Bounceback scheme is triggered
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=Inf; % Purosely assigned infinity, later on will be replaced with correct value
            else
                ;
            end
        end
    elseif NInt==2031 % avg bounceback of f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2032 % local average Rho + linear extrapolated f_neq
        ;
    elseif NInt==2033 % linear extrapolated Rho + linear extrapolated f_neq
        ;
    elseif NInt==2098 % Given x-velocity and Rho, y velociy is extrapolated
        U_ex=node_star(ND,U,0);
        f_bc_eq(:,1)=eqm_h(V,[U_bc(1);U_ex(2)],Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NInt==2099 % Given U and Rho on boundary + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
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
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on Density Inlet boundary is not valid or available!');
    end
elseif abs(ND{2})==30 % Velocity Outlet
    if ND{2}<0
        error('Immersed boundary cannot be velocity outlet!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    if NOut==3099
%         f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc=f_bc_eq+f_bc_neq;
        
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
%         F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
%             F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on Velocity Outlet boundary is not valid or available!');
    end
elseif abs(ND{2})==31 % Pressure Outlet
    if ND{2}<0
        error('Immersed boundary cannot be pressure outlet!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    % f = f_eq + f_neq
    if NOut==3099 % Given U and Rho on boundary + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3100 % local average U + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3101 % linear extrapolated U + local average f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,2*U_i-U_ii,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3102 % local average U + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3103 % linear extrapolated U + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,2*U_i-U_ii,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + bounce-back f_neq
    elseif NOut==3105 % f = f_eq + f_neq. Second-order mapping for the unknow part 
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_U=zeros(num_star_cell,2);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_U(i,:)=U(:,star_cell_index(i))';
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,([1,ND{3}']*(CoeX*F_U))',Rho_bc,T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3130 % local average Rho + local average f_neq
        B=1; % Bounceback scheme is triggered
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=Inf; % Purosely assigned infinity, later on will be replaced with correct value
            else
                ;
            end
        end
    elseif NOut==3131 % avg bounceback of f_neq
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NOut==3199 % v fixed to zero & extrapolated u + local average f_neq
%         U_ex=node_star(ND,U,0);
%         f_bc_eq(:,1)=eqm_h(V,[U_ex(1);0],Rho_bc,T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc=f_bc_eq+f_bc_neq;

        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
%         F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
%             F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on Density Outlet boundary is not valid or available!');
    end
elseif abs(ND{2})==4 % Stationary Wall
    bc=ND{21};
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    % f = f_eq
    if NStw==400 % fixed density
        f_bc(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
    elseif NStw==401 % local average Rho
        f_bc(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
    elseif NStw==402 % linear extrapolated Rho
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        % f = f_eq + f_neq
    elseif NStw==410 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc_neq(:,1)=node_star(ND,f-feq,0);
        % BB &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%         if V==V1
%             UB=ND{19};
%         elseif V==V2
%             UB=ND{20};
%         else
%             error('The lattice used is not the lattice used for determing upwind cells');
%         end
%         f_bc_neq_temp=f_bc_neq;
%         for l=1:q
%             if UB(l)~=l
%                 f_bc_neq_temp(l,1)=f_bc_neq(UB(l),1); % replaced with correct number
%             end
%         end
%         f_bc_neq=f_bc_neq_temp;
        % BB &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        if fdcy==0 % No decay
            f_bc=f_bc_eq+f_bc_neq;
        elseif fdcy==1 % with decay
            f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
        else
            error('The flag for decaying of f_neq is invalid!');
        end
    elseif NStw==411 % linear extrapolated Rho + local average f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==412 % local average Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==413 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with fixed density
    elseif NStw==414 % assembled Rho from f_node + linear extrapolated f_neq, This one proves not working. see NODE v02.06.2019
        cell_stencil=ND{23};
        Coef_stencil=ND{25};
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
%         f_bc_eq(:,1)=eqm_h(V,U_bc,RHO(1,cell_stencil(1,1)),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=f(:,cell_stencil(1,1))-feq(:,cell_stencil(1,1));
%         f_bc_eq(:,1)=eqm_h(V,U_bc,Coef_stencil(1,1)*RHO(1,cell_stencil(1,1))+Coef_stencil(1,2)*RHO(1,cell_stencil(1,2)),T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc_neq(:,1)=Coef_stencil(1,1)*(f(:,cell_stencil(1,1))-feq(:,cell_stencil(1,1)))+Coef_stencil(1,2)*(f(:,cell_stencil(1,2))-feq(:,cell_stencil(1,2)));
        f_bc=f_bc_eq+f_bc_neq;
%         f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         [Rho_nd_neq,U_nd_neq,T_nd_neq]=macro_h(f_bc_neq,V,2,fd);
%         if Rho_nd_neq<0
%             disp('Rho_nd_neq<0');
%         end
%         f_bc_neq_relx(:,1)=eqm_h(V,U_bc,Rho_nd_neq,T_bc,q,wh,Rho_r,fd);
%         f_bc=f_bc_eq+f_bc_neq_relx;
        
%         [Rho_nd,U_nd,T_nd]=macro_h(f_node(:,ND{1}),V,2,fd);
%         f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_nd,T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc=f_node(:,ND{1})-1/tau*f_bc_neq;
    elseif NStw==415 % f = f_eq + f_neq. Second-order mapping for the unknow part
        
        %% The following scheme that is commented does not work due to instability . See NODE v04.11.2019
%         cell_center=CELL{ND{23}};
%         cell_neighbor=cell_center{38};
%         CoeA=cell_center{40};
%         
%         F_RHO=[(RHO(:,cell_neighbor(1))-RHO(:,ND{23}))';(RHO(:,cell_neighbor(2))-RHO(:,ND{23}))';(RHO(:,cell_neighbor(3))-RHO(:,ND{23}))'];
%         f_bc_eq(:,1)=eqm_h(V,U_bc,RHO(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_RHO))',T_bc,q,wh,Rho_r,fd);
%         
%         F_f=[(f(:,cell_neighbor(1))-f(:,ND{23}))';(f(:,cell_neighbor(2))-f(:,ND{23}))';(f(:,cell_neighbor(3))-f(:,ND{23}))'];
%         F_feq=[(feq(:,cell_neighbor(1))-feq(:,ND{23}))';(feq(:,cell_neighbor(2))-feq(:,ND{23}))';(feq(:,cell_neighbor(3))-feq(:,ND{23}))'];
%         f_bc_neq(:,1)=(f(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_f))')-(feq(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_feq))');
        
        
        
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,([1,ND{3}']*(CoeX*F_RHO))',T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==420 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==421 % linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % Bounceback rules
    elseif NStw==430 % f = f_eq + bounce-back only for unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        %         if abs(f_bc_neq(2,1))>abs(f_bc_neq(2,1))
        %             f_bc_neq(2,1)=f_bc_neq(4,1);
        %         else
        %             f_bc_neq(4,1)=f_bc_neq(2,1);
        %         end
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==431 % bounceback of all f_neq (switch opposing components)
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==432  % f = f_eq + avg bouncebacked all f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==433 % f = f_eq + avg bouncebacked only unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_neq=f_bc_neq_temp;
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq_temp(UB(l),1); % Bounceback of opposing component
                f_bc_neq(UB(l),1)=f_bc_neq_temp(l,1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==434 % f = f_eq + bouncebacked only the normal unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        % f_bc_neq=f_neq_extp/2
    elseif NStw==440
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq/2;
    elseif NStw==441 % with bounceback
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=(node_star(ND,f,0)-node_star(ND,feq,0))/2;
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        % With collision
    elseif NStw==450 % no bb
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
    elseif NStw==451 % bb of only the normal unknown pdf, derivative of 434
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+(1-1/tau)*f_bc_neq;
        % BB + counter-slip
    elseif NStw==460  % BB of only the normal unknown pdf + counter slip of f_neq, derivatibe of 434
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        [Rho_c,U_c]=macro_h(f_bc_neq,V,Rho_r,fd);
        if Rho_c<0
            error('The counter-slip density is negative!');
        end
        f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
        f_bc_neq=f_bc_neq-f_c;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NStw==461  % BB of only the normal unknown pdf + counter slip of f_bc with f_c determined by maxiwilliam, derivatibe of 434
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            U_c=U_bb-U_bc;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc-f_c;
        else
            Rho_c=Rho_bc-Rho_bb;
            U_c=U_bc-U_bb;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc+f_c;
        end
    elseif NStw==462  % BB of only the normal unknown pdf + counter slip of f_bc with f_c=Rho_c and only added to the component that is bouncebacked, derivatibe of 434
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)-Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        else
            Rho_c=Rho_bc-Rho_bb;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)+Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        end
    elseif NStw==463  % BB of all unknown pdf + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
    elseif NStw==464  % BB of all unknown f + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=node_star(ND,f,0);
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % BB + counter-slip by R.S. Miller + collision
    elseif NStw==470  % collision before counter-slip
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % Evolution at nodes
    elseif NStw==480 %
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=f_node(:,ND{1})-dt/tau*(f_node(:,ND{1})-eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd));
    elseif NStw==499 % both velocity and density are prescibed
        Rho_bc=bc(1,1);
%         f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
%         f_bc=f_bc_eq+f_bc_neq;

        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
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
    % f = f_eq
    if NMow==500 % fixed density
        f_bc(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
    elseif NMow==501 % local average Rho
        f_bc(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
    elseif NMow==502 % linear extrapolated Rho
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        % f = f_eq + f_neq
    elseif NMow==510 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        if fdcy==0 % No decay
            f_bc=f_bc_eq+f_bc_neq;
        elseif fdcy==1 % with decay
            f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
        else
            error('The flag for decaying of f_neq is invalid!');
        end
    elseif NMow==511 % linear extrapolated Rho + local average f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,~]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==512 % local average Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==513 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with fixed density
    elseif NMow==515 % f = f_eq + f_neq. Second-order mapping for the unknow part
        %% The following scheme that is commented does not work due to instability . See NODE v04.11.2019
%         cell_center=CELL{ND{23}};
%         cell_neighbor=cell_center{38};
%         CoeA=cell_center{40};
%         
%         F_RHO=[(RHO(:,cell_neighbor(1))-RHO(:,ND{23}))';(RHO(:,cell_neighbor(2))-RHO(:,ND{23}))';(RHO(:,cell_neighbor(3))-RHO(:,ND{23}))'];
%         f_bc_eq(:,1)=eqm_h(V,U_bc,RHO(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_RHO))',T_bc,q,wh,Rho_r,fd);
%         
%         F_f=[(f(:,cell_neighbor(1))-f(:,ND{23}))';(f(:,cell_neighbor(2))-f(:,ND{23}))';(f(:,cell_neighbor(3))-f(:,ND{23}))'];
%         F_feq=[(feq(:,cell_neighbor(1))-feq(:,ND{23}))';(feq(:,cell_neighbor(2))-feq(:,ND{23}))';(feq(:,cell_neighbor(3))-feq(:,ND{23}))'];
%         f_bc_neq(:,1)=(f(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_f))')-(feq(:,ND{23})+((ND{3}-cell_center{5})'*(CoeA*F_feq))');
        
        
        
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_RHO=zeros(num_star_cell,1);
        F_f=zeros(num_star_cell,qh);
        F_feq=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_RHO(i,:)=RHO(1,star_cell_index(i));
            F_f(i,:)=f(:,star_cell_index(i))';
            F_feq(i,:)=feq(:,star_cell_index(i))';
        end
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,([1,ND{3}']*(CoeX*F_RHO))',T_bc,q,wh,Rho_r,fd);

        f_bc_neq(:,1)=([1,ND{3}']*(CoeX*F_f))'-([1,ND{3}']*(CoeX*F_feq))';

        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==520 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==521 % linear extrapolated f_neq
        NS=ND{16};
        NDI=NODE{NS(3)};
        NDII=NODE{NS(4)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + bounce-back f_neq
    elseif NMow==530 % f = f_eq + bounce-back only for unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
            if abs(f_bc_neq(2,1))>abs(f_bc_neq(2,1))
                f_bc_neq(2,1)=f_bc_neq(4,1);
            else
                f_bc_neq(4,1)=f_bc_neq(2,1);
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==531 % bounceback of all f_neq (switch opposing components)
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==532  % f = f_eq + avg bouncebacked f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==533 % f = f_eq + avg bouncebacked only unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_neq=f_bc_neq_temp;
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq_temp(UB(l),1); % Bounceback of opposing component
                f_bc_neq(UB(l),1)=f_bc_neq_temp(l,1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==534 % f = f_eq + bouncebacked only the normal unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        % f_bc_neq=f_neq_extp/2
    elseif NMow==540
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq/2;
    elseif NMow==541 % with bounceback
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=(node_star(ND,f,0)-node_star(ND,feq,0))/2;
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        % With collision
    elseif NMow==550 % no bb
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
    elseif NMow==551 % bb of only the normal unknown pdf, derivative of 534
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+(1-1/tau)*f_bc_neq;
        % BB + counter-slip
    elseif NMow==560  % BB of only the normal unknown pdf + counter slip of f_neq with f_c determined by maxiwilliam, derivatibe of 534
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        [Rho_c,U_c]=macro_h(f_bc_neq,V,Rho_r,fd);
        if Rho_c<0
            error('The counter-slip density is negative!');
        end
        f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
        f_bc_neq=f_bc_neq-f_c;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NMow==561  % BB of only the normal unknown pdf + counter slip of f_bc with f_c determined by maxiwilliam, derivatibe of 534
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            U_c=U_bb-U_bc;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc-f_c;
        else
            Rho_c=Rho_bc-Rho_bb;
            U_c=U_bc-U_bb;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc+f_c;
        end
    elseif NMow==562  % BB of only the normal unknown pdf + counter slip of f_bc with f_c=Rho_c and only added to the component that is bouncebacked, derivatibe of 534
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        Ncu=NDU{2};
        Ncd=NDD{2};
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)-Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        else
            Rho_c=Rho_bc-Rho_bb;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)+Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        end
    elseif NMow==563  % BB of all unknown f_neq + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
    elseif NMow==564  % BB of all unknown f + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=node_star(ND,f,0);
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % BB + counter-slip by R.S. Miller + collision
    elseif NMow==570  % collision before counter-slip
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % Evolution at nodes
    elseif NMow==580 %
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=f_node(:,ND{1})-dt/tau*(f_node(:,ND{1})-eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd));
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==6 % Well Developed
    if ND{2}<0
        error('Immersed boundary cannot be Well Developed!');
    end
    % Zero gradient, f = f_upwind
    if NWdp==600
        NP=ND{16};
        f_bc=f_node(:,NP(3));
        % Linear gradient, f = 2*f_upwind - f_up_upwind
    elseif NWdp==610
        NP=ND{16};
        f_bc=2*f_node(:,NP(3))-f_node(:,NP(4));
        % Linear gradient, f = node_star
    elseif NWdp==615
        num_star_cell=ND{4};
        star_cell_index=ND{5};
        CoeX=ND{23};
        F_f=zeros(num_star_cell,qh);
        for i=1:num_star_cell
            F_f(i,:)=f(:,star_cell_index(i))';
        end

        f_bc=([1,ND{3}']*(CoeX*F_f))';
    elseif NWdp==620
        f_bc=node_star(ND,f,0);
        % f = f_eq + f_neq, linear extrapolation of Rho, U
    elseif NWdp==630
        NP=ND{16};
        U_e=2*U_node(:,NP(3))-U_node(:,NP(4));
        RHO_e=2*RHO_node(:,NP(3))-RHO_node(:,NP(4));
        f_bc_eq(:,1)=eqm_h(V,U_e,RHO_e,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NWdp==640 % This is place holder to match up fort the same code number in pdf_bc_t
        f_bc=node_star(ND,f,0);
    else
        error('Input scheme code for nodal pdf on Well Developed boundary is not valid or available!');
    end
elseif abs(ND{2})==8 % Inflow or Outflow
    if ND{2}<0
        error('Immersed boundary cannot be Inflow or Outflow!');
    end
    if NIof==800
        f_bc=node_star(ND,f,0);
    else
        error('Input scheme code for nodal pdf on Inflow/Outflow boundary is not valid or available!');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Corner nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif abs(ND{2})==70 % Periodic + Periodic
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Periodic!');
    end
    if NC70==7000
        f_bc=node_star(ND,f,1);
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
    % f = f_eq
    if NC71==7100 % local average density
        f_bc(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
    elseif NC71==7101 % Weighted linear extrapolated Rho in two vertical directions
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        Rho_p=0.5*(Rho_L+Rho_R);
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC71!');
        end
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        Rho_np=2*Rho_i-Rho_ii;
        f_bc(:,1)=eqm_h(V,U_bc,wc*Rho_p+(1-wc)*Rho_np,q,wh,Rho_r,fd);
    elseif NC71==7102 % linear extrapolated Rho along diagonal direction
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        % f = f_eq + f_neq with weighted two-side extrapolation
    elseif NC71==7110 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7111 % weighted linear extrapolated Rho + local average f_neq
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        Rho_p=0.5*(Rho_L+Rho_R);
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC71!');
        end
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        Rho_np=2*Rho_i-Rho_ii;
        f_bc_eq(:,1)=eqm_h(V,U_bc,wc*Rho_p+(1-wc)*Rho_np,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7112 % local average Rho + weighted linear extrapolated f_neq
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_neq_p(:,1)=0.5*((f_node(:,NDL{1})-f_eq_l(:,1))+(f_node(:,NDR{1})-f_eq_r(:,1)));
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC71!');
        end
        f_eq_i(:,1)=node_star(NDI,feq,1);
        f_eq_ii(:,1)=node_star(NDII,feq,1);
        f_bc_neq_np(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_p(:,1)+(1-wc)*f_bc_neq_np(:,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7113 % weighted linear extrapolated Rho + weighted linear extrapolated f_neq
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        [Rho_L,U_L]=macro_h(f_node(:,NDL{1}),V,Rho_r,fd);
        [Rho_R,U_R]=macro_h(f_node(:,NDR{1}),V,Rho_r,fd);
        Rho_p=0.5*(Rho_L+Rho_R);
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_neq_p(:,1)=0.5*((f_node(:,NDL{1})-f_eq_l(:,1))+(f_node(:,NDR{1})-f_eq_r(:,1)));
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC71!');
        end
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        Rho_np=2*Rho_i-Rho_ii;
        f_eq_i(:,1)=node_star(NDI,feq,1);
        f_eq_ii(:,1)=node_star(NDII,feq,1);
        f_bc_neq_np(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,wc*Rho_p+(1-wc)*Rho_np,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_p(:,1)+(1-wc)*f_bc_neq_np(:,1);
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with extrapolation along diagonal
    elseif NC71==7120 % local average Rho + local average f_neq SAME AS NC71=7110
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7121 % linear extrapolated Rho + local average f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7122 % local average Rho + linear extrapolated f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7123 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with fixed density
    elseif NC71==7130 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7131 % weighted linear extrapolated f_neq
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_neq_p(:,1)=0.5*((f_node(:,NDL{1})-f_eq_l(:,1))+(f_node(:,NDR{1})-f_eq_r(:,1)));
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC71!');
        end
        f_eq_i(:,1)=node_star(NDI,feq,1);
        f_eq_ii(:,1)=node_star(NDII,feq,1);
        f_bc_neq_np(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_p(:,1)+(1-wc)*f_bc_neq_np(:,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7132 % f_neq with extrapolation along diagonal
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_total
    elseif NC71==7140 % copy of exterior neighbor
        NS=ND{16};
        f_bc(:,1)=f_node(:,NS(1));
    elseif NPdc==7141 % local average
        f_bc=node_star(ND,f,1);
        % Bounceback rule
    elseif NC71==7150 % f = f_eq + bounce-back only for unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7152  % f = f_eq + avg bouncebacked f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7153 % f = f_eq + avg bouncebacked only unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc_neq=f_bc_neq_temp;
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq_temp(UB(l),1); % Bounceback of opposing component
                f_bc_neq(UB(l),1)=f_bc_neq_temp(l,1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7154 % f = f_eq + bouncebacked only normal unknown f_neq
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        if NDU{2}==1 && NDD{2}~=1
            Ncu=ND{2};
            Ncd=NDD{2};
        elseif NDU{2}~=1 && NDD{2}==1
            Ncu=NDU{2};
            Ncd=ND{2};
        else
            error('The current corner node should have been the joint of periodic and non-periodic BC!');
        end
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        % with collision
    elseif NC71==7160 % Collision without BB
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+(1-1/tau)*f_bc_neq;
    elseif NC71==7161 % BB of only the normal unknown pdf, derivative of 7154
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        if NDU{2}==1 && NDD{2}~=1
            Ncu=ND{2};
            Ncd=NDD{2};
        elseif NDU{2}~=1 && NDD{2}==1
            Ncu=NDU{2};
            Ncd=ND{2};
        else
            error('The current corner node should have been the joint of periodic and non-periodic BC!');
        end
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,1),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+(1-1/tau)*f_bc_neq;
        % BB + counter-slip
    elseif NC71==7170  % BB of only the normal unknown pdf + counter slip of f_neq, derivatibe of 7154
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        if NDU{2}==1 && NDD{2}~=1
            Ncu=ND{2};
            Ncd=NDD{2};
        elseif NDU{2}~=1 && NDD{2}==1
            Ncu=NDU{2};
            Ncd=ND{2};
        else
            error('The current corner node should have been the joint of periodic and non-periodic BC!');
        end
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        [Rho_c,U_c]=macro_h(f_bc_neq,V,Rho_r,fd);
        if Rho_c<0
            error('The counter-slip density is negative!');
        end
        f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
        f_bc_neq=f_bc_neq-f_c;
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC71==7171  % BB of only the normal unknown pdf + counter slip of f_bc with f_c determined by maxiwilliam, derivatibe of 7154
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        if NDU{2}==1 && NDD{2}~=1
            Ncu=ND{2};
            Ncd=NDD{2};
        elseif NDU{2}~=1 && NDD{2}==1
            Ncu=NDU{2};
            Ncd=ND{2};
        else
            error('The current corner node should have been the joint of periodic and non-periodic BC!');
        end
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            U_c=U_bb-U_bc;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc-f_c;
        else
            Rho_c=Rho_bc-Rho_bb;
            U_c=U_bc-U_bb;
            f_c=eqm_h(V,U_c(:,1),Rho_c(1,:),qh,wh,Rho_r,fd);
            f_bc=f_bc+f_c;
        end
    elseif NC71==7172  % BB of only the normal unknown pdf + counter slip of f_bc with f_c=Rho_c and only added to the pdf that is bouncebacked, derivatibe of 7154
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Nng=ND{16};
        NDU=NODE{Nng(1)};
        NDD=NODE{Nng(2)};
        if NDU{2}==1 && NDD{2}~=1
            Ncu=ND{2};
            Ncd=NDD{2};
        elseif NDU{2}~=1 && NDD{2}==1
            Ncu=NDU{2};
            Ncd=ND{2};
        else
            error('The current corner node should have been the joint of periodic and non-periodic BC!');
        end
        L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
        n_x=-(Ncd(2)-Ncu(2))/L;
        n_y=(Ncd(1)-Ncu(1))/L;
        
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        a=0;
        for l=1:q
            if UB(l)~=l
                if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                    f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
                    a=1;
                end
            else
                ;
            end
        end
        if a==0
            error('Bounce-back is not applied correctly!');
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        if Rho_bb>Rho_bc
            Rho_c=Rho_bb-Rho_bc;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)-Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        else
            Rho_c=Rho_bc-Rho_bb;
            a=0;
            for l=1:q
                if UB(l)~=l
                    if [n_x,n_y]*V(:,l)<0 && single(e+n_x*V(2,l)-n_y*V(1,l))==single(e)
                        f_bc(l,1)=f_bc(l,1)+Rho_c; % Bounceback of opposing component
                        a=1;
                    end
                else
                    ;
                end
            end
            if a==0
                error('Bounce-back is not applied correctly!');
            end
        end
    elseif NC71==7173  % BB of all unknown pdf + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,1);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==72 % Periodic + Density or Density + Periodic
    bc=ND{21};
    Rho_bc=bc(1,1);
    T_bc=bc(4,1);
    if ND{2}<0
        error('Immersed boundary cannot be Periodic + Density or Density + Periodic!');
    end
    if NC72==7200
        ;
    elseif NC72==7210 % Local Rho + extrapolated U
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,1),Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,1)-node_star(ND,feq,1);
        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on periodic + density boundary is not valid or available!');
    end
elseif abs(ND{2})==73 % Velocity + Velocity
    bc=ND{21};
    Rho_bc=bc(1,1);
    U_bc=bc(2:3,1);
    T_bc=bc(4,1);
    % f = f_eq
    if NC73==7300 % local average density
        f_bc(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
    elseif NC73==7301 % Weighted linear extrapolated Rho in two vertical directions
        NS=ND{16};
        NDI_U=NODE{NS(1)};
        NDI_D=NODE{NS(2)};
        NSU=NDI_U{16};
        NSD=NDI_D{16};
        NDII_U=NODE{NSU(1)};
        NDII_D=NODE{NSD(2)};
        [Rho_i_u,U_i_u]=macro_h(f_node(:,NDI_U{1}),V,Rho_r,fd);
        [Rho_ii_u,U_ii_u]=macro_h(f_node(:,NDII_U{1}),V,Rho_r,fd);
        Rho_u=2*Rho_i_u-Rho_ii_u;
        [Rho_i_d,U_i_d]=macro_h(f_node(:,NDI_D{1}),V,Rho_r,fd);
        [Rho_ii_d,U_ii_d]=macro_h(f_node(:,NDII_D{1}),V,Rho_r,fd);
        Rho_d=2*Rho_i_d-Rho_ii_d;
        f_bc(:,1)=eqm_h(V,U_bc,wc*Rho_u+(1-wc)*Rho_d,q,wh,Rho_r,fd);
    elseif NC73==7302 % linear extrapolated Rho along diagonal direction
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        % f = f_eq + f_neq with weighted two-side extrapolation
    elseif NC73==7310 % local average Rho + local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        if fdcy==0 % No decay
            f_bc=f_bc_eq+f_bc_neq;
        elseif fdcy==1 % with decay
            f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
        else
            error('The flag for decaying of f_neq is invalid!');
        end
    elseif NC73==7311 % weighted linear extrapolated Rho + local average f_neq
        NS=ND{16};
        NDI_U=NODE{NS(1)};
        NDI_D=NODE{NS(2)};
        NSU=NDI_U{16};
        NSD=NDI_D{16};
        NDII_U=NODE{NSU(1)};
        NDII_D=NODE{NSD(2)};
        [Rho_i_u,U_i_u]=macro_h(f_node(:,NDI_U{1}),V,Rho_r,fd);
        [Rho_ii_u,U_ii_u]=macro_h(f_node(:,NDII_U{1}),V,Rho_r,fd);
        Rho_u=2*Rho_i_u-Rho_ii_u;
        [Rho_i_d,U_i_d]=macro_h(f_node(:,NDI_D{1}),V,Rho_r,fd);
        [Rho_ii_d,U_ii_d]=macro_h(f_node(:,NDII_D{1}),V,Rho_r,fd);
        Rho_d=2*Rho_i_d-Rho_ii_d;
        f_bc_eq(:,1)=eqm_h(V,U_bc,wc*Rho_u+(1-wc)*Rho_d,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7312 % local average Rho + weighted linear extrapolated f_neq
        NS=ND{16};
        NDI_U=NODE{NS(1)};
        NDI_D=NODE{NS(2)};
        NSU=NDI_U{16};
        NSD=NDI_D{16};
        NDII_U=NODE{NSU(1)};
        NDII_D=NODE{NSD(2)};
        f_eq_i_u(:,1)=node_star(NDI_U,feq,0);
        f_eq_ii_u(:,1)=node_star(NDII_U,feq,0);
        f_bc_neq_u(:,1)=2*(f_node(:,NDI_U{1})-f_eq_i_u(:,1))-(f_node(:,NDII_U{1})-f_eq_ii_u(:,1));
        f_eq_i_d(:,1)=node_star(NDI_D,feq,0);
        f_eq_ii_d(:,1)=node_star(NDII_D,feq,0);
        f_bc_neq_d(:,1)=2*(f_node(:,NDI_D{1})-f_eq_i_d(:,1))-(f_node(:,NDII_D{1})-f_eq_ii_d(:,1));
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_u(:,1)+(1-wc)*f_bc_neq_d(:,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7313 % weighted linear extrapolated Rho + weighted linear extrapolated f_neq
        NS=ND{16};
        NDI_U=NODE{NS(1)};
        NDI_D=NODE{NS(2)};
        NSU=NDI_U{16};
        NSD=NDI_D{16};
        NDII_U=NODE{NSU(1)};
        NDII_D=NODE{NSD(2)};
        [Rho_i_u,U_i_u]=macro_h(f_node(:,NDI_U{1}),V,Rho_r,fd);
        [Rho_ii_u,U_ii_u]=macro_h(f_node(:,NDII_U{1}),V,Rho_r,fd);
        Rho_u=2*Rho_i_u-Rho_ii_u;
        [Rho_i_d,U_i_d]=macro_h(f_node(:,NDI_D{1}),V,Rho_r,fd);
        [Rho_ii_d,U_ii_d]=macro_h(f_node(:,NDII_D{1}),V,Rho_r,fd);
        Rho_d=2*Rho_i_d-Rho_ii_d;
        f_eq_i_u(:,1)=node_star(NDI_U,feq,0);
        f_eq_ii_u(:,1)=node_star(NDII_U,feq,0);
        f_bc_neq_u(:,1)=2*(f_node(:,NDI_U{1})-f_eq_i_u(:,1))-(f_node(:,NDII_U{1})-f_eq_ii_u(:,1));
        f_eq_i_d(:,1)=node_star(NDI_D,feq,0);
        f_eq_ii_d(:,1)=node_star(NDII_D,feq,0);
        f_bc_neq_d(:,1)=2*(f_node(:,NDI_D{1})-f_eq_i_d(:,1))-(f_node(:,NDII_D{1})-f_eq_ii_d(:,1));
        f_bc_eq(:,1)=eqm_h(V,U_bc,wc*Rho_u+(1-wc)*Rho_d,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_u(:,1)+(1-wc)*f_bc_neq_d(:,1);
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with extrapolation along diagonal
    elseif NC73==7314 % assembled Rho from f_node + collision on bc
        cell_stencil=ND{23};
        Coef_stencil=ND{25};
        f_bc_eq(:,1)=eqm_h(V,U_bc,Coef_stencil(1,1)*RHO(1,cell_stencil(1,1))+Coef_stencil(1,2)*RHO(1,cell_stencil(1,2)),T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=Coef_stencil(1,1)*(f(:,cell_stencil(1,1))-feq(:,cell_stencil(1,1)))+Coef_stencil(1,2)*(f(:,cell_stencil(1,2))-feq(:,cell_stencil(1,2)));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7320 % local average Rho + local average f_neq SAME AS NC73=7110
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7321 % linear extrapolated Rho + local average f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7322 % local average Rho + linear extrapolated f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7323 % linear extrapolated Rho + linear extrapolated f_neq
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        [Rho_i,U_i]=macro_h(f_node(:,NDI{1}),V,Rho_r,fd);
        [Rho_ii,U_ii]=macro_h(f_node(:,NDII{1}),V,Rho_r,fd);
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,2*Rho_i-Rho_ii,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_eq + f_neq with fixed density
    elseif NC73==7330 % local average f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7331 % weighted linear extrapolated f_neq
        NS1=ND{16};
        NDL=NODE{NS1(1)};
        NDR=NODE{NS1(3)};
        f_eq_l(:,1)=node_star(NDL,feq,0);
        f_eq_r(:,1)=node_star(NDR,feq,0);
        f_bc_neq_p(:,1)=0.5*((f_node(:,NDL{1})-f_eq_l(:,1))+(f_node(:,NDR{1})-f_eq_r(:,1)));
        NS2=ND{16};
        NUM=NODE{NS2(1)};
        NDM=NODE{NS2(2)};
        if NUM{2}==1 && NDM{2}~=1
            NDI=NUM;
            NS=NDI{16};
            NDII=NODE{NS(1)};
        elseif NUM{2}~=1 && NDM{2}==1
            NDI=NDM;
            NS=NDI{16};
            NDII=NODE{NS(2)};
        else
            error('The current corner node is not NC73!');
        end
        f_eq_i(:,1)=node_star(NDI,feq,1);
        f_eq_ii(:,1)=node_star(NDII,feq,1);
        f_bc_neq_np(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_p(:,1)+(1-wc)*f_bc_neq_np(:,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7332 % f_neq with extrapolation along diagonal
        NS=ND{9};
        NDI=NODE{NS(1)};
        NDII=NODE{NS(2)};
        f_eq_i(:,1)=node_star(NDI,feq,0);
        f_eq_ii(:,1)=node_star(NDII,feq,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=2*(f_node(:,NDI{1})-f_eq_i(:,1))-(f_node(:,NDII{1})-f_eq_ii(:,1));
        f_bc=f_bc_eq+f_bc_neq;
        % f = f_total
    elseif NC73==7340
        ;
    elseif NC73==7341 % local average
        f_bc=node_star(ND,f,0);
        % Transient model, solve LBM on zero CV at boundary node
    elseif NC73==7350
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc=(1-dt/tau)*f_node(:,ND{1})+dt/tau*f_bc_eq(:,1);
    elseif NC73==7351
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_r,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq(:,1)+(1-dt/tau)*f_bc_neq(:,1);
        % Bounceback rule
    elseif NC73==7360  % f = f_eq + bouncebacked f_neq
        %         if V==V1
        %             UB=ND{19};
        %         elseif V==V2
        %             UB=ND{20};
        %         else
        %             error('The lattice used is not the lattice used for determing upwind cells');
        %         end
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC73==7361 % bounceback of total pdf
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_temp=f_bc_eq+f_bc_neq;
        %% Temperoraly for D2Q9
        f_bc(1,1)=f_bc_temp(1,1);
        f_bc(2,1)=f_bc_temp(4,1);
        f_bc(4,1)=f_bc_temp(2,1);
        f_bc(3,1)=f_bc_temp(5,1);
        f_bc(5,1)=f_bc_temp(3,1);
        f_bc(6,1)=f_bc_temp(8,1);
        f_bc(8,1)=f_bc_temp(6,1);
        f_bc(7,1)=f_bc_temp(9,1);
        f_bc(9,1)=f_bc_temp(7,1);
    elseif NC73==7362  % f = f_eq + avg bouncebacked f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
        % Bounceback + counter slip
    elseif NC73==7370  % BB of all unknown pdf + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
    elseif NC73==7371  % BB of all unknown pdf + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bb; % Only difference from 7370
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
    elseif NC73==7372  % BB of all unknown f + counter slip of f_bc by R.S. Miller
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=node_star(ND,f,0);
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % with collision
    elseif NC73==7380 % no bb
        f_bc_eq(:,1)=eqm_h(V,U_bc,node_star(ND,RHO,0),q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
    elseif NC73==7381 % bb + counter-slip + collision before counter-slip
        if V==V1
            UB=ND{19};
        elseif V==V2
            UB=ND{20};
        else
            error('The lattice used is not the lattice used for determing upwind cells');
        end
        Rho_bc=node_star(ND,RHO,0);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        for l=1:q
            if UB(l)~=l
                f_bc_neq(l,1)=f_bc_neq(UB(l),1); % Bounceback of opposing component
            else
                ;
            end
        end
        f_bc=f_bc_eq+(1-dt/tau)*f_bc_neq;
        [Rho_bb,U_bb]=macro_h(f_bc,V,Rho_r,fd);
        Mm=(U_bb-U_bc)*Rho_bc;
        alpha=0.5;
        for l=1:q
            if UB(l)~=l
                f_bc(l,1)=f_bc(l,1)-(Mm')*V(:,l)*alpha; % Redistribute mass
            else
                ;
            end
        end
        % Evolution at nodes
    elseif NC73==7390 %
        Rho_bc=node_star(ND,RHO,0);
        f_bc(:,1)=f_node(:,ND{1})-dt/tau*(f_node(:,ND{1})-eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd));
        % Special edition
    elseif NC73==7399 % both velocity and density are prescibed
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on Moving Wall boundary is not valid or available!');
    end
elseif abs(ND{2})==74 % Density + Density
    if ND{2}<0
        error('Immersed boundary cannot be Density + Density!');
    end
    bc=ND{21};
    Rho_bc=bc(1,1);
    T_bc=bc(4,1);
    if NC74==7400 % local average density
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Input scheme code for nodal pdf on Density + Density corner node is not valid or available!');
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
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC75==7501 % weighted linear extrapolated f_neq
        NS=ND{16};
        NDI_U=NODE{NS(1)};
        NDI_D=NODE{NS(2)};
        NSU=NDI_U{16};
        NSD=NDI_D{16};
        NDII_U=NODE{NSU(1)};
        NDII_D=NODE{NSD(2)};
        f_eq_i_u(:,1)=node_star(NDI_U,feq,0);
        f_eq_ii_u(:,1)=node_star(NDII_U,feq,0);
        f_bc_neq_u(:,1)=2*(f_node(:,NDI_U{1})-f_eq_i_u(:,1))-(f_node(:,NDII_U{1})-f_eq_ii_u(:,1));
        f_eq_i_d(:,1)=node_star(NDI_D,feq,0);
        f_eq_ii_d(:,1)=node_star(NDII_D,feq,0);
        f_bc_neq_d(:,1)=2*(f_node(:,NDI_D{1})-f_eq_i_d(:,1))-(f_node(:,NDII_D{1})-f_eq_ii_d(:,1));
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=wc*f_bc_neq_u(:,1)+(1-wc)*f_bc_neq_d(:,1);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC75==7502 % local average Rho and U
%         f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),node_star(ND,RHO,0),T_bc,q,wh,Rho_r,fd);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,T_bc,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    elseif NC75==7510  % f = f_eq + avg bouncebacked f_neq
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,q,wh,Rho_r,fd);
        f_bc_neq_temp(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        %% Temperoraly for D2Q9
        f_bc_neq(1,1)=f_bc_neq_temp(1,1);
        f_bc_neq(2,1)=f_bc_neq_temp(4,1);
        f_bc_neq(4,1)=f_bc_neq_temp(2,1);
        f_bc_neq(3,1)=f_bc_neq_temp(5,1);
        f_bc_neq(5,1)=f_bc_neq_temp(3,1);
        f_bc_neq(6,1)=f_bc_neq_temp(8,1);
        f_bc_neq(8,1)=f_bc_neq_temp(6,1);
        f_bc_neq(7,1)=f_bc_neq_temp(9,1);
        f_bc_neq(9,1)=f_bc_neq_temp(7,1);
        f_bc_neq=(f_bc_neq+f_bc_neq_temp)/2;
        f_bc=f_bc_eq+f_bc_neq;
    end
elseif abs(ND{2})==76 % Well Developed + others or others + Well Developed
    if ND{2}<0
        error('Immersed boundary cannot be Well Developed + others or others + Well Developed!');
    end
    % Zero gradient, f = f_upwind
    if NC76==7600
        NS=ND{18};
        NUM=NODE{NS(1)};
        NDM=NODE{NS(2)};
        if NUM{2}~=6
            f_bc=f_node(:,NS(1));
        elseif NDM{2}~=6
            f_bc=f_node(:,NS(2));
        else
            error('Check the flag of corner node!')
        end
        % Linear gradient, f = 2*f_upwind - f_upwind_upwind
    elseif NC76==7601
        NS=ND{18};
        NUM=NODE{NS(1)};
        NDM=NODE{NS(2)};
        if NUM{2}~=6
            f_bc=f_node(:,NS(2));
        elseif NDM{2}~=6
            f_bc=f_node(:,NS(1));
        else
            error('Check the flag of corner node!')
        end
    elseif NC76==7603
        NS=ND{18};
        NUM=NODE{NS(1)};
        NDM=NODE{NS(2)};
        if NUM{2}~=6
            ND_center=NUM;
        elseif NDM{2}~=6
            ND_center=NDM;
        else
            error('Check the flag of corner node!')
        end
        NSS=ND_center{18};
        NUUM=NODE{NSS(1)};
        NDDM=NODE{NSS(2)};
        if NUUM{2}~=76
            ND_left=NUUM;
        elseif NDDM{2}~=76
            ND_left=NDDM;
        else
            error('Check the flag of corner node!')
        end
        f_bc=2*f_node(:,ND_center{1})-f_node(:,ND_left{1});
    elseif NC76==7610
        NS=ND{18};
        NUM=NODE{NS(1)};
        NDM=NODE{NS(2)};
        NSU=NUM{16};
        NSD=NDM{16};
        if NUM{2}~=6
            f_bc=2*f_node(:,NS(1))-f_node(:,NSU(1));
        elseif NDM{2}~=6
            f_bc=2*f_node(:,NS(2))-f_node(:,NSD(2));
        else
            error('Check the flag of corner node!')
        end
        % Zero gradient, f = node_star
    elseif NC76==7620
        f_bc=node_star(ND,f,0);
        % f = f_eq + f_neq, linear extrapolation of Rho, U
    elseif NC76==7621
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),node_star(ND,RHO,0),1,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_neq(:,1)=node_star(ND,f-feq,0);
        f_bc=f_bc_eq+f_bc_neq;
%         f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
    elseif NC76==7622
        bc=ND{21};
        Rho_bc=bc(1,1);
        f_bc_eq(:,1)=eqm_h(V,node_star(ND,U,0),Rho_bc,1,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_neq(:,1)=node_star(ND,f-feq,0);
        f_bc=f_bc_eq+f_bc_neq;
%         f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
    elseif NC76==7623 % for both density and velocity are known
        bc=ND{21};
        Rho_bc=bc(1,1);
        U_bc=bc(2:3,1);
        f_bc_eq(:,1)=eqm_h(V,U_bc,Rho_bc,1,q,wh,Rho_r,fd);
%         f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc_neq(:,1)=node_star(ND,f-feq,0);
        f_bc=f_bc_eq+f_bc_neq;
%         f_bc=f_bc_eq+f_bc_neq*exp(-dt/tau);
    elseif NC76==7630 %% This is abondoned
        NS=ND{18};
        NUM=NODE{NS(1)};
        NDM=NODE{NS(2)};
        NSU=NUM{18};
        NSD=NDM{18};
        if NUM{2}~=6
            U_e=2*U_node(:,NS(1))-U_node(:,NSU(1));
            RHO_e=2*RHO_node(:,NS(1))-RHO_node(:,NSU(1));
        elseif NDM{2}~=6
            U_e=2*U_node(:,NS(2))-U_node(:,NSD(2));
            RHO_e=2*RHO_node(:,NS(2))-RHO_node(:,NSD(2));
        else
            error('Check the flag of corner node!')
        end
        f_bc_eq(:,1)=eqm_h(V,U_e,RHO_e,q,wh,Rho_r,fd);
        f_bc_neq(:,1)=node_star(ND,f,0)-node_star(ND,feq,0);
        f_bc=f_bc_eq+f_bc_neq;
    else
        error('Check the flag of corner node!');
    end
elseif abs(ND{2})==77 % Well Developed + Well Developed
    if NC77==7700
        f_bc=node_star(ND,f,0);
    else
        error('Check the flag of corner node!')
    end
elseif abs(ND{2})==78 % Periodic + Inflow/Outflow
    if NC78==7800
        f_bc=node_star(ND,f,1);
    else
        error('Check the flag of corner node!')
    end
elseif abs(ND{2})==79 % Velocity or Density + Inflow/Outflow
    if NC79==7900
        f_bc=node_star(ND,f,0);
    else
        error('Check the flag of corner node!')
    end
elseif abs(ND{2})==80 % Inflow/Outflow + Inflow/Outflow
    if NC80==8000
        f_bc=node_star(ND,f,0);
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