%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%PDF Initialization%%%%%%%%%%%%%%%%%%%%%%%%%
switch FF
    case 0
        %@@@@@@@@@@@@@@@@@@@@@@ 0. Uniform flow @@@@@@@@@@@@@@@@@@@@@@@@@@
        %%%Initial Condition of cell centriods
        if strcmp(top,'Moving Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % Couette flow, analytic solution is known
            %%%Initial Condition of cell centriods
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
                %%% Calculate PDF
                f_old(1:qh,r)=eqm_h(V,[u_x_ana;u_y_ana],Rho_r,t_ana,qh,wh,Rho_r,FD);
            end;
            f_old_old=f_old;
            %%%%Initial Condition of nodes
            for r=1:N;
                Node=NODE{r};
                if Node{2}~=0 && Node{2}~=1
                    bc_node=Node{21};
                    f_nd_old(1:qh,r)=eqm_h(V,bc_node(2:3),Rho_r,bc_node(4),qh,wh,Rho_r,FD);
                else
                    f_nd_old(:,r)=node_star(Node,f_old,0);
                end
            end
            f_nd=f_nd_old;
            f_nd_eq=f_nd;
        elseif strcmp(left,'Velocity Inlet')==1 % General Poiselle flow
            f_old=zeros(qh,M);
            f_nd_old=zeros(qh,N);
            if strcmp(top,'Stationary Wall')==1 || strcmp(bottom,'Stationary Wall')==1 % Parabolic velocity profile
                %%%Initial Condition of cell centriods
                for r=1:M;
                    Cell=CELL{r};
                    coor=Cell{5};
                    %%% Calculate PDF
                    f_old(1:qh,r)=eqm_h(V,[sym_para(Y1,(Y2-Y1),U_in(2,4),coor(1,1));0],Rho_r,1,qh,wh,Rho_r,FD);
                end;
                f_old_old=f_old;
                %%%%Initial Condition of nodes
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
            else % Same velocity profile as inlet
                %%%Initial Condition of cell centriods
                for r=1:M;
                    Cell=CELL{r};
                    coor=Cell{5};
                    %%% Calculate PDF
                    if F_u_in(1,4)==2
                        f_old(1:qh,r)=eqm_h(V,[sym_para(Y1,(Y2-Y1),U_in(2,4),coor(1,1));0],Rho_r,1,qh,wh,Rho_r,FD);
                    else
                        f_old(1:qh,r)=eqm_h(V,[U_in(2,4);0],Rho_r,1,qh,wh,Rho_r,FD);
                    end
                end;
                f_old_old=f_old;
                %%%%Initial Condition of nodes
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
            end
        else
            %%%Initial Condition of cell centriods
            for r=1:M;
                f_old(1:qh,r)=eqm_h(V,[0;0],Rho_r,1,qh,wh,Rho_r,FD);
            end
            f_old_old=f_old;
            %%%%Initial Condition of nodes
            for r=1:N;
                f_nd_old(1:qh,r)=eqm_h(V,[0;0],Rho_r,1,qh,wh,Rho_r,FD);
                f_nd=f_nd_old;
                f_nd_eq=f_nd;
            end
        end
    case 1
        %@@@@@@@@@@@@@@@@@@@@ 1. Taylor vortex flow @@@@@@@@@@@@@@@@@@@@@@@
        k1=2*pi/(X2-X1);
        k2=2*pi/(Y2-Y1);
        U_0=0.1;
        %%%Initial Condition of cell centriods
        for r=1:M
            P=CELL{r};
            Centroid=P{5};
            f_old(1:qh,r)=eqm_h(V,[-U_0*cos(k1*Centroid(1,1))*sin(k2*Centroid(2,1));U_0*(k1/k2)*sin(k1*Centroid(1,1))*cos(k2*Centroid(2,1))],Rho_r,1,qh,wh,Rho_r,FD);
        end
        %%%%Initial Condition of nodes
        for r=1:N
            f_nd_old(1:qh,r)=eqm_h(V,[-U_0*cos(k1*X(r))*sin(k2*Y(r));U_0*(k1/k2)*sin(k1*X(r))*cos(k2*Y(r))],Rho_r,1,qh,wh,Rho_r,FD);
            f_nd=f_nd_old;
        end
        f_nd_eq=f_nd;
    case 2
        %@@@@@@@@@@@@@@@@@@@@ 2. Perturbation flow @@@@@@@@@@@@@@@@@@@@@@@
        %%%Initial Condition of cell centriods
        for r=1:M;
            f_old(1:qh,r)=eqm_h(V,[0;0],Rho_r,1,qh,wh,Rho_r,FD);
        end;
        %%%%Initial Condition of nodes
        for r=1:N;
            f_nd_old(1:qh,r)=eqm_h(V,[0;0],Rho_r,1,qh,wh,Rho_r,FD);
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
[Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
%%%%Initialize density & velocity at cnodes
[Rho_nd,U_nd,T_nd]=macro_h(f_nd_old,V,Rho_r,FD);
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
%%