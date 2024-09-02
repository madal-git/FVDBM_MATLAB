%%%% This is function to initialize thermal PDF for Passive-Scalar mdoel          

if FTH==1 && qh~=37
    switch FF
        case 0
            %@@@@@@@@@@@@@@@@@@@@@@ 0. Uniform flow @@@@@@@@@@@@@@@@@@@@@@@@@@
            %%%Initial Condition of cell centriods
            if strcmp(top,'Moving Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % Couette flow, analytic solution is known
                %%%Initial Condition of cell centriods
                g_old=zeros(qt,M);
                for r=1:M;
                    Cell=CELL{r};
                    coor=Cell{5};
                    %%% Analytic velocity solution
                    u_x_ana=U_m(1,1)/(Y2-Y1)*coor(2,1);
                    u_y_ana=0;
                    %%% Analytic Temperature solution
                    Pr=1; % For BGK and single set of PDF for both hydrodynamics and thermodynamics
                    Cp=2; % Cp=(D/2)+1
                    Ec=U_m(1,1)^2/Cp/(T_bc(1)-T_bc(3));
                    t_ana=((coor(2,1)/(Y2-Y1))+Pr*Ec/2*((coor(2,1)/(Y2-Y1))*(1-(coor(2,1)/(Y2-Y1)))))*(T_bc(1)-T_bc(3))+T_bc(3);
                    %%% Calculate PDF for hydrodynamics
                    g_old(:,r)=eqm_t(V2,[u_x_ana;u_y_ana],Rho_r,t_ana,qt,wt,Rho_r,FD);
                end;
                g_old_old=g_old;
                g=g_old;
                g_eq=g;
                g_eq_old=g_eq;
                %%%%Initial Condition of nodes
                g_nd_old=zeros(qt,N);
                for r=1:N;
                    Node=NODE{r};
                    if (Node{2}~=0 && Node{2}~=1) && Node{2}~=6
                        bc_node=Node{21};
                        g_nd_old(:,r)=eqm_t(V2,bc_node(2:3),bc_node(1),bc_node(4),qt,wt,Rho_r,FD);
                    else
                        g_nd_old(:,r)=node_star(Node,g_old,0);
                    end
                end
                g_nd=g_nd_old;
                g_nd_eq=g_nd;
            else
                error('Other type of flows are temporarily not available for thermal model');
            end
    end
end


T=macro_t(g_old,V2,Rho);
%%%%Initialize density & velocity at cnodes
T_nd=macro_t(g_nd_old,V2,Rho_nd);
%%%%%%%%%%Density & Velocity Initialization based on PDF%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Initialization of Monitered Variables%%%%%%%%%%%%%%%
%########################## 1. Global ############################
%%%% Initial value of monitor density and velocity and residual
%%%% PDF residual initialization
R_t(1)=1;
%%%% Initial value of monitor Integral flux
FLX_c_t=zeros(qt,1);