%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%PDF Initialization%%%%%%%%%%%%%%%%%%%%%%%%%
switch FF
    case 0
        %@@@@@@@@@@@@@@@@@@@@@@ 0. Uniform flow @@@@@@@@@@@@@@@@@@@@@@@@@@
        %%%Initial Condition of cell centriods
        if strcmp(top,'Moving Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % Couette flow, analytic solution is known
            %%%Initial Condition of cell centriods
            for r=1:M
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
                f_old(1:qh,r)=eqm_h(V,[u_x_ana;u_y_ana],Rho_ref,t_ana,qh,wh,Rho_ref,FD);
            end
            f_old_old=f_old;
            %%%%Initial Condition of nodes
            for r=1:N
                Node=NODE{r};
                if (Node{2}~=0 && Node{2}~=1) && Node{2}~=6
                    bc_node=Node{21};
                    f_nd_old(1:qh,r)=eqm_h(V,bc_node(2:3),bc_node(1),bc_node(4),qh,wh,Rho_ref,FD);
                else
                    f_nd_old(:,r)=node_star(Node,f_old,0);
                end
            end
            f_nd=f_nd_old;
            f_nd_eq=f_nd;
        elseif (strcmp(top,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1) && (strcmp(left,'Pressure Inlet')==1 && strcmp(right,'Pressure Outlet')==1) % Pressure driven flow Poiselle flow
                Rho_l=Rho_in(4);
                Rho_r=Rho_out(2);
                a=(Rho_l-Rho_r)/(X1-X2);
                b=Rho_l-a*X1;
                if qh==9
                    dpdx=(Rho_l-Rho_r)/(X2-X1)/3;
                    Mew=Tau/3*Rho_ref;
                    U_ana_max_out=dpdx/Mew/2*((Y2-Y1)/2)^2;
                else
                    error('Other lattice is not available1');
                end
                %%%Initial Condition of cell centriods
                for r=1:M;
                    Cell=CELL{r};
                    coor=Cell{5};
                    %%% Calculate PDF
                    f_old(1:qh,r)=eqm_h(V,[U_ana_max_out*(1-((coor(2)-(Y2-Y1)/2)/((Y2-Y1)/2))^2);0],a*coor(1)+b,1,qh,wh,Rho_ref,FD);
                end
                f_old_old=f_old;
                %%%%Initial Condition of nodes
                f_nd_old=zeros(qh,N);
                for r=1:N;
                    Node=NODE{r};
                    if Node{2}~=0
                        ;
                    elseif Node{2}==1
                        f_nd_old(:,r)=node_star(Node,f_old,1);
                    else
                        f_nd_old(:,r)=node_star(Node,f_old,0);
                    end
                end
                f_nd=f_nd_old;
                f_nd_eq=f_nd;
        elseif (strcmp(top,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1) && (strcmp(left,'Velocity Inlet')==1) % Velocity driven Poiselle flow
                f_old=zeros(qh,M);
                f_nd_old=zeros(qh,N);
                    %%%Initial Condition of cell centriods
                    for r=1:M
                        Cell=CELL{r};
                        coor=Cell{5};
                        %%% Calculate PDF
                        f_old(1:qh,r)=eqm_h(V,[sym_para(Y1,(Y2-Y1),U_in(1,4),coor(2,1));0],Rho_ref,1,qh,wh,Rho_ref,FD);
                    end
                    f_old_old=f_old;
                    %%%%Initial Condition of nodes
                    for r=1:N
                        Node=NODE{r};
                        if Node{2}~=0
                            ;
                        elseif Node{2}==1
                            f_nd_old(:,r)=node_star(Node,f_old,1);
                        else
                            f_nd_old(:,r)=node_star(Node,f_old,0);
                        end
                    end
                    f_nd=f_nd_old;
                    f_nd_eq=f_nd;
        elseif strcmp(left,'Pressure Inlet')==1 && strcmp(right,'Pressure Outlet')==1  % General Pressure driven flow
                Rho_l=Rho_in(4);
                Rho_r=Rho_out(2);
                a=(Rho_l-Rho_r)/(X1-X2);
                b=Rho_l-a*X1;
                %%%Initial Condition of cell centriods
                for r=1:M
                    Cell=CELL{r};
                    coor=Cell{5};
                    %%% Calculate PDF
                    f_old(1:qh,r)=eqm_h(V,[0;0],a*coor(1)+b,1,qh,wh,Rho_ref,FD);
                end
                f_old_old=f_old;
                %%%%Initial Condition of nodes
                f_nd_old=zeros(qh,N);
                for r=1:N
                    Node=NODE{r};
                    if Node{2}~=0
                        ;
                    elseif Node{2}==1
                        f_nd_old(:,r)=node_star(Node,f_old,1);
                    else
                        f_nd_old(:,r)=node_star(Node,f_old,0);
                    end
                end
                f_nd=f_nd_old;
                f_nd_eq=f_nd;
        elseif strcmp(top,'Pressure Outlet')==1 && strcmp(bottom,'Pressure Outlet')==1 % Stratified flow
%             %%% Fill density data for the right boundary
%             for r=1:N
%                 if (r>N_I+N_L-1 && r<N_I+N_L-1+N_H-1) || (r>N-N_H+1 && r<N) % The right boundary and the left boundary
%                     Node=NODE{r};
%                     coor=Node{3};
%                     bc_node=Node{21};
%                     bc_node(1)=(coor(2,1)-Y1)*(Rho_out(1,3)-Rho_in(1,1))/(Y1-Y2)+Rho_out(1,3);
%                     Node{21}=bc_node;
%                     NODE{r}=Node;
%                 end
%             end
            %%%Initial Condition of cell centriods
            for r=1:M
                Cell=CELL{r};
                coor=Cell{5};
                rho=(coor(2,1)-Y1)*(Rho_out(1,3)-Rho_out(1,1))/(Y1-Y2)+Rho_out(1,3);
                %%% Calculate PDF for hydrodynamics
                f_old(1:qh,r)=eqm_h(V,U_in(:,4),rho,T_ref,qh,wh,Rho_ref,FD);
            end
            f_old_old=f_old;
            %%%%Initial Condition of nodes
            for r=1:N
                Node=NODE{r};
%                 if (Node{2}~=0 && Node{2}~=1) && (Node{2}~=6 && Node{2}~=76)
%                     bc_node=Node{21};
%                     f_nd_old(1:qh,r)=eqm_h(V,bc_node(2:3),bc_node(1),bc_node(4),qh,wh,Rho_ref,FD);
%                 else
                    f_nd_old(:,r)=node_star(Node,f_old,0);
%                 end
            end
            f_nd=f_nd_old;
            f_nd_eq=f_nd;
        elseif (strcmp(top,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1) && (strcmp(left,'Stationary Wall')==1 && strcmp(right,'Stationary Wall')==1) % Stratified flow in a box
%             %%% Fill density data for the right boundary
%             for r=1:N
%                 if (r>N_I+N_L-1 && r<N_I+N_L-1+N_H-1) || (r>N-N_H+1 && r<N) % The right boundary and the left boundary
%                     Node=NODE{r};
%                     coor=Node{3};
%                     bc_node=Node{21};
%                     bc_node(1)=(coor(2,1)-Y1)*(Rho_out(1,3)-Rho_in(1,1))/(Y1-Y2)+Rho_out(1,3);
%                     Node{21}=bc_node;
%                     NODE{r}=Node;
%                 end
%             end
            %%%Initial Condition of cell centriods
            for r=1:M
                Cell=CELL{r};
                coor=Cell{5};
                rho=(coor(2,1)-Y1)*(Rho_out(1,3)-Rho_out(1,1))/(Y1-Y2)+Rho_out(1,3);
                %%% Calculate PDF for hydrodynamics
                f_old(1:qh,r)=eqm_h(V,[0;0],rho,T_ref,qh,wh,Rho_ref,FD);
            end
            f_old_old=f_old;
            %%%%Initial Condition of nodes
            for r=1:N
                Node=NODE{r};
                if Node{2}~=0
                    bc_node=Node{21};
                    f_nd_old(1:qh,r)=eqm_h(V,bc_node(2:3),bc_node(1),bc_node(4),qh,wh,Rho_ref,FD);
                else
                    f_nd_old(:,r)=node_star(Node,f_old,0);
                end
            end
            f_nd=f_nd_old;
            f_nd_eq=f_nd;
        else
            %%%Initial Condition of cell centriods
            for r=1:M
%                 f_old(1:qh,r)=eqm_h(V,[0;0],Rho_ref,1,qh,wh,Rho_ref,FD);
%                 % Zero velocity
                f_old(1:qh,r)=eqm_h(V,(rand(2,1)-0.5)/100,Rho_ref,1,qh,wh,Rho_ref,FD); % Pertubation
            end
            f_old_old=f_old;
            f_eq=f_old;
            %%%%Initial Condition of nodes
            for r=1:N
                f_nd_old(1:qh,r)=eqm_h(V,[0;0],Rho_ref,1,qh,wh,Rho_ref,FD);
                f_nd=f_nd_old;
                f_nd_eq=f_nd;
            end

%             for r=1:M
%                 Cell=CELL{r};
%                 coor=Cell{5};
%                 %%% Calculate PDF
%                 %                     if F_u_in(1,4)==2
%                 %                         f_old(1:qh,r)=eqm_h(V,[sym_para(Y1,(Y2-Y1),U_in(1,4),coor(2,1));0],Rho_ref,1,qh,wh,Rho_ref,FD);
%                 %                     else
%                 %                         f_old(1:qh,r)=eqm_h(V,[U_in(1,4);0],Rho_ref,1,qh,wh,Rho_ref,FD);
%                 %                     end
%                 %                         f_old(1:qh,r)=eqm_h(V,[0;0],Rho_ref,1,qh,wh,Rho_ref,FD); % Zero velocity
%                 f_old(1:qh,r)=eqm_h(V,(rand(2,1)-0.5)/100,Rho_ref,1,qh,wh,Rho_ref,FD); % Pertubation
%             end
%             f_old_old=f_old;
%             f_eq=f_old;
%             %%%%Initial Condition of nodes
%             for r=1:N
%                 Node=NODE{r};
%                 %                     if Node{2}~=0
%                 %                         ;
%                 %                     elseif Node{2}==1
%                 %                         f_nd_old(:,r)=node_star(Node,f_old,1);
%                 %                     else
%                 f_nd_old(:,r)=node_star(Node,f_old,0);
%                 %                     end
%             end
%             f_nd=f_nd_old;
%             f_nd_eq=f_nd;
        end
    case 1
        %@@@@@@@@@@@@@@@@@@@@ 1. Taylor vortex flow @@@@@@@@@@@@@@@@@@@@@@@
        k1=2*pi/(X2-X1);
        k2=2*pi/(Y2-Y1);
        if qh==7
            c_s=(1/sqrt(4));
        elseif qh==9
            c_s=(1/sqrt(3));
        elseif qh==13
            c_s=(1/sqrt(2));
        else
            error('Wrong lattice!');
        end
        U_0=0.01*c_s;
        %%%Initial Condition of cell centriods
        for r=1:M
            P=CELL{r};
            Centroid=P{5};
            f_old(1:qh,r)=eqm_h(V,[-U_0*cos(k1*Centroid(1,1))*sin(k2*Centroid(2,1));U_0*(k1/k2)*sin(k1*Centroid(1,1))*cos(k2*Centroid(2,1))],Rho_ref+(-U_0^2/4*(cos(2*k1*Centroid(1,1))+k1^2/k2^2*cos(2*k2*Centroid(2,1))))/c_s^2,1,qh,wh,Rho_ref,FD);
        end
        %%%%Initial Condition of nodes
        for r=1:N
            f_nd_old(1:qh,r)=eqm_h(V,[-U_0*cos(k1*X(r))*sin(k2*Y(r));U_0*(k1/k2)*sin(k1*X(r))*cos(k2*Y(r))],Rho_ref+(-U_0^2/4*(cos(2*k1*X(r))+k1^2/k2^2*cos(2*k2*Y(r))))/c_s^2,1,qh,wh,Rho_ref,FD);
            f_nd=f_nd_old;
        end
        f_nd_eq=f_nd;
        t_c=log(2)/(Mew1*(k1^2+k2^2))
    case 2
        %@@@@@@@@@@@@@@@@@@@@ 2. Perturbation flow @@@@@@@@@@@@@@@@@@@@@@@
        %%%Initial Condition of cell centriods
        for r=1:M;
            f_old(1:qh,r)=eqm_h(V,[0;0],Rho_ref,1,qh,wh,Rho_ref,FD);
        end;
        %%%%Initial Condition of nodes
        for r=1:N;
            f_nd_old(1:qh,r)=eqm_h(V,[0;0],Rho_ref,1,qh,wh,Rho_ref,FD);
            f_nd=f_nd_old;
        end;
        %%%%%%%%% Introduce perturbation to symmetric four triangles at the
        %%%%%%%%% center of the square, if not needed, comment this section.
        PERT=0.1; %%%% Pertubation of density
        XP=(X1+X2)/2;
        YP=(Y1+Y2)/2;
        if FM==0
            spc=h;
        elseif FM==1
            spc=(dx+dy)/2;
        else
            error('Wrong flag for mesh type!');
        end
        for l=1:N
            ND=NODE{l};
            if dis(ND{3},[XP;YP])<spc/10
                break;
            end
        end
        RAND=ND{5};
        for i=1:ND{4}
            %%%%%% Mass pertubation
            f_old(1:qh,RAND(i))=f_old(1:qh,RAND(i))+PERT/(qh-1)*ones(qh,1);
            %%%%%% Momentum pertubation
            %     f_old(2,RAND(i))=f_old(2,RAND(i))+PERT/(qh-1);
            %     f_old(4,RAND(i))=f_old(4,RAND(i))-PERT/(qh-1);
        end
end
%%%%%%%%%%%%%%%%%%%%%%%PDF Initialization%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Density & Velocity Initialization based on PDF%%%%%%%%%%
%%%%Initialize density & velocity at cell centriods
[Rho,U,T]=macro_h(f_old,V,Rho_ref,FD);
%%%%Initialize density & velocity at cnodes
[Rho_nd,U_nd,T_nd]=macro_h(f_nd_old,V,Rho_ref,FD);
%%%%%%%%%%Density & Velocity Initialization based on PDF%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Initialization of Monitered Variables%%%%%%%%%%%%%%%
%########################## 1. Global ############################
%%%% Initial value of monitor density and velocity and residual
RHO(1)=0;
for r=1:M
    RHO(1)=RHO(1)+Rho(1,r);
end
RHO(1)=RHO(1)/M;
    
UR(:,1)=[0;0];
for r=1:M
    UR(:,1)=UR(:,1)+U(:,r);
end
UR(:,1)=UR(:,1)/M;

TR(1)=0;
for r=1:M
    TR(1)=TR(1)+T(1,r);
end
TR(1)=TR(1)/M;
%%%% PDF residual initialization
R(1)=1;
%%%% Initial value of monitor Integral flux
FLX_c=zeros(qh,1);
%%%% Initial value of monitor Integral collision
FCOL_c=zeros(qh,1);
%%%%%%%%%%% Finding the checking points
%########################## 2. Local ############################%
% if FM==0
%     for r=1:N
%         if single(X(r))==single((X2-X1)/2) && single(Y(r))==single((Y2-Y1)/2);
%             Tar=r;
%             break;
%         end
%     end
%     spc=[1];
%     a=length(spc);
%     NDX=zeros(2,a);
%     NDY=zeros(2,a);
%     for r=1:N
%         for l=1:a
%             if single(X(r))==single(X(Tar)+spc(l)*h) && single(Y(r))==single(Y(Tar))
%                 NDX(1,l)=r;
%             elseif single(X(r))==single(X(Tar)-spc(l)*h) && single(Y(r))==single(Y(Tar))
%                 NDX(2,l)=r;
%             elseif single(X(r))==single(X(Tar)) && single(Y(r))==single(Y(Tar)+spc(l)*h)
%                 NDY(1,l)=r;
%             elseif single(X(r))==single(X(Tar)) && single(Y(r))==single(Y(Tar)-spc(l)*h)
%                 NDY(2,l)=r;
%             end
%         end
%     end
%     %%%%%%%%%% Initialize the checking variables
%     for l=1:a
%         RHO_NDX(1+2*(l-1):2*l,1)=[Rho_nd(NDX(1,l));Rho_nd(NDX(2,l))];
%         RHO_NDY(1+2*(l-1):2*l,1)=[Rho_nd(NDY(1,l));Rho_nd(NDY(2,l))];
%         U_NDX(1+4*(l-1):4*l,1)=[U_nd(:,NDX(1,l));U_nd(:,NDX(2,l))];
%         U_NDY(1+4*(l-1):4*l,1)=[U_nd(:,NDY(1,l));U_nd(:,NDY(2,l))];
%     end
%     
%     if FF==1 % Taylor vortex flow
%         for l=1:a
%             U_TV_NDX(1+4*(l-1):4*l,1)=[U_nd(:,NDX(1,l));U_nd(:,NDX(2,l))];
%             U_TV_NDY(1+4*(l-1):4*l,1)=[U_nd(:,NDY(1,l));U_nd(:,NDY(2,l))];
%         end
%     end
% end
%%%%%%%%%%%%%%Initialization of Monitered Variables%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thermal_constant=8.314;
%% Initialization for thermal model
if FTH==1
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
            if qh==qt
                g_old(:,r)=t_ana*f_old(:,r);
            else
                g_old(:,r)=eqm_t(V2,[u_x_ana;u_y_ana],Rho_ref,t_ana,qt,wt,Rho_ref,FD);
            end
        end;
    elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % General horizontal Poiseulle flow
        g_old=zeros(qt,M);
        % For Poiseulle flow
%         for r=1:M
%             Cell=CELL{r};
%             coor=Cell{5};
%             T(1,r)=sin(pi*coor(2,1)/(Y2-Y1))+T_ref;
%         end
        % For natural convection between two parallel plates
        for r=1:M
            Cell=CELL{r};
            coor=Cell{5};
            T(1,r)=T_bc(3)+(T_bc(3)-T_bc(1))/(Y1-Y2)*(coor(2,1)-Y1)+((rand(1,1)-0.5)*T_ref)/100;
        end
        for r=1:M
%             if qh==qt
% % %                 g_old(:,r)=T_bc(1)*f_old(:,r);
% %             Cell=CELL{r};
% %             coor=Cell{5};
% %             T_ana=T_bc(3)+(T_bc(1)-T_bc(3))/(Y2-Y1)*coor(2,1);
% %             g_old(:,r)=T_ana*f_old(:,r);
%                 g_old(:,r)=Thermal_constant*T_bc(1)*f_old(:,r);
%             else
                g_old(:,r)=eqm_t(V2,U(:,r),Rho_ref,T(1,r),qt,wt,Rho_ref,FD);
%             end
        end
    elseif strcmp(top,'Periodic')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Periodic')==1 && strcmp(left,'Periodic')==1 % periodic flow
        g_old=zeros(qt,M);
        if FCD==1 %% Pure conduction
            %%% Calculate PDF for hydrodynamics
            for r=1:M
                Rho(1,r)=Rho_ref;
                U(:,r)=[0;0];
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T_ref,qh,wh,Rho_ref,FD);
                f_old(:,r)=f_eq(:,r);
            end
            T_max=1;
            for r=1:M
                Cell=CELL{r};
                coor=Cell{5};
                T(1,r)=T_max*exp(-((coor(1,1)-(X2-X1)/2)^2+(coor(2,1)-(Y2-Y1)/2)^2)/((X2-X1)/4+(Y2-Y1)/4)^2);
            end
            for r=1:M
                g_old(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_ref,FD);
            end
        else
            error('Other type of flows are temporarily not available for thermal model');
        end
    elseif strcmp(top,'Fully Developed')==1 && strcmp(right,'Fully Developed')==1 && strcmp(bottom,'Fully Developed')==1 && strcmp(left,'Fully Developed')==1 % periodic flow
        g_old=zeros(qt,M);
        if FCD==1 %% Pure conduction
            %%% Calculate PDF for hydrodynamics
            for r=1:M
                Rho(1,r)=Rho_ref;
                U(:,r)=[0;0];
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T_ref,qh,wh,Rho_ref,FD);
                f_old(:,r)=f_eq(:,r);
            end
            T_max=1;
            for r=1:M
                Cell=CELL{r};
                coor=Cell{5};
                T(1,r)=T_max*exp(-((coor(1,1)-(X2-X1)/2)^2+(coor(2,1)-(Y2-Y1)/2)^2)/((X2-X1)/4+(Y2-Y1)/4)^2);
            end
            for r=1:M
                g_old(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_ref,FD);
            end
        else
            error('Other type of flows are temporarily not available for thermal model');
        end
    elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1 % periodic flow
        g_old=zeros(qt,M);
        if FCD==1 %% Pure conduction
            %%% Calculate PDF for hydrodynamics
            for r=1:M
                Rho(1,r)=Rho_ref;
                U(:,r)=[0;0];
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T_ref,qh,wh,Rho_ref,FD);
                f_old(:,r)=f_eq(:,r);
            end
            T_max=1;
            for r=1:M
                Cell=CELL{r};
                coor=Cell{5};
                T(1,r)=T_max*exp(-((coor(1,1)-(X2-X1)/2)^2+(coor(2,1)-(Y2-Y1)/2)^2)/((X2-X1)/4+(Y2-Y1)/4)^2);
            end
            for r=1:M
                g_old(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_ref,FD);
            end
        else
            error('Other type of flows are temporarily not available for thermal model');
        end
    elseif N_I_N>0 % convection around inner boundary
        g_old=zeros(qt,M);
        for r=1:M
%             Cell=CELL{r};
%             coor=Cell{5};
            if qh==qt
% %                 g_old(:,r)=T_bc(1)*f_old(:,r);
%             Cell=CELL{r};
%             coor=Cell{5};
%             T_ana=T_bc(3)+(T_bc(1)-T_bc(3))/(Y2-Y1)*coor(2,1);
%             g_old(:,r)=T_ana*f_old(:,r);
                g_old(:,r)=T_ref*f_old(:,r);
            else
                g_old(:,r)=eqm_t(V2,[0;0],Rho_ref,T_bc(4),qt,wt,Rho_ref,FD);
            end
        end
    else
        error('Other type of flows are temporarily not available for thermal model');
    end
    g_old_old=g_old;
    g=g_old;
    g_eq=g;
    g_eq_old=g_eq;
    %% Initial Condition of nodes
    g_nd_old=zeros(qt,N);
    for r=1:N
        Node=NODE{r};
        if (Node{2}~=0 && Node{2}~=1) && (Node{2}~=6 && Node{2}~=77)
            if qh==qt
                bc_nd=Node{21};
                g_nd_old(:,r)=Thermal_constant*bc_nd(4,1)*f_nd_old(:,r);
            else
                bc_node=Node{21};
                g_nd_old(:,r)=eqm_t(V2,bc_node(2:3),bc_node(1),bc_node(4),qt,wt,Rho_ref,FD);
            end
        else
            g_nd_old(:,r)=node_star(Node,g_old,0);
        end
    end
    g_nd=g_nd_old;
    g_nd_eq=g_nd;
    %% Macro for T
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
end


[Rho,U,T]=macro_h(f_old,V,Rho_ref,FD);
%%%% Calculate equilibrium pdf at cell centroids
for r=1:M
    f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_ref,FD);
end
for l=1:N
    NP=NODE{l};
    if NP{2}==0 % Interior nodes
        % if W~=1
        f_nd(:,l)=node_star(NP,f_old,0);
        % end
    else % nodes on boundary
        f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_ref,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80);
    end
end
[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_ref,FD);
if FTH==1 && FTHM~=2
    T=macro_t(g_old,V2,Rho);
    %%%% Calculate equilibrium pdf at cell centroids
%     if qt==qh
%         for r=1:M
%             g_eq(:,r)=T(1,r)*f_eq(:,r);
%         end
%     else
        for r=1:M
            g_eq(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_ref,FD);
        end
%     end
    for l=1:N
        NP=NODE{l};
        if NP{2}==0 % Interior nodes
            g_nd(:,l)=node_star(NP,g_old,0);
        else % nodes on boundary
            g_nd(:,l)=pdf_bc_t(N_L,N_H,N_I,Tau_t,dt,NODE,NP,f_old,f_eq,f_nd,g_old,g_eq,g_nd,U_nd,U,Rho_nd,Rho,V2,V1,V2,Rho_ref,Rho_in,Rho_out,qh,qt,wc,wh,wt,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80);
        end
    end
    T_nd=macro_t(g_nd,V2,Rho_nd);
end
%%% Calculate the velocity megnitude in cell centroid
U_M=zeros(1,M);
for r=1:M
    P=CELL{r};
    U_M(r)=sqrt(U(1,r)^2+U(2,r)^2);
end
%%% Adding the boundary values in the matrices for display
XXX=zeros(1,M);
YYY=zeros(1,M);
%%%% coordinates
for r=1:M
    P=CELL{r};
    Centroid=P{5};
    XXX(r)=Centroid(1,1);
    YYY(r)=Centroid(2,1);
end
for r=1:N-N_I
    XXX(M+r)=X(r+N_I);
    YYY(M+r)=Y(r+N_I);
end
for r=1:N_I_N
    XXX(M+N-N_I+r)=X(r);
    YYY(M+N-N_I+r)=Y(r);
end
U_plt=U;
Rho_plt=Rho;
T_plt=T;
for r=1:N-N_I
    U_plt(:,M+r)=U_nd(:,r+N_I);
    U_M(M+r)=sqrt(U_nd(1,r+N_I)^2+U_nd(2,r+N_I)^2);
    Rho_plt(1,M+r)=Rho_nd(1,r+N_I);
    T_plt(1,M+r)=T_nd(1,r+N_I);
end
for r=1:N_I_N
    U_plt(:,M+N-N_I+r)=U_nd(:,r);
    U_M(M+N-N_I+r)=sqrt(U_nd(1,r)^2+U_nd(2,r)^2);
    Rho_plt(1,M+N-N_I+r)=Rho_nd(1,r);
    T_plt(1,M+N-N_I+r)=T_nd(1,r);
end