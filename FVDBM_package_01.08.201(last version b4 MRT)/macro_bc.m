function NODE=macro_bc(X1,X2,Y1,Y2,N_I,N_H,N_L,NODE,Rho_r,Rho_in,F_rho_in,Rho_out,F_rho_out,U_r,U_m,U_in,F_u_in,U_out,F_u_out,T_bc,TG_bc,F_T,F_TG,FCO)
% function NODE=macro_bc(NODE,Rho_r,Rho_in,F_rho_in,Rho_out,F_rho_out,U_m,U_in,F_u_in,U_out,F_u_out,T_in,T_out,TG_in,TG_out)
% writes all macroscopic variables into NODE data structure
% NODE - The node data structure
% Rho_in - The density at inlet
% Rho_out - The density at outlet
% U_in - The velocity at inlet
% U_out - The velocity at outlet
% U_m - The velocity at moving wall
% T_bc - The temperature at all boundary walls
% TG_bc - The temperature gradient at all boundary walls

N=length(NODE);

%% Hydro
% Check validity of U_m
U_top=U_m(:,1);
U_right=U_m(:,2);
U_bottom=U_m(:,3);
U_left=U_m(:,4);
if (U_top(2,1)~=0 || U_bottom(2,1)~=0) || (U_left(1,1)~=0 || U_right(1,1)~=0)
    error('Check U_m for moving wall velocity!');
end
% Check positivity of density
for i=length(Rho_r)
    if Rho_r(i)<=0
        error('Density has to be positive!');
    end
end
for i=length(Rho_in)
    if Rho_in(i)<=0
        error('Density has to be positive!');
    end
end
for i=length(Rho_out)
    if Rho_out(i)<=0
        error('Density has to be positive!');
    end
end
% Hydro_bc_corner is used to store hydro bc parameters in corner nodes when
% flag of corner node is 73 or 74
Lc=length(U_r(:,1))+1; % Total 4 positions, [Rho;Ux;Uy;T]
Hydro_bc_corner=zeros(2*Lc,N);
% Pre-fill NP{21} with zero vector
for l=1:N
    NP=NODE{l};
    NP{21}=zeros(Lc,1);
    NODE{l}=NP;
end
% Write hydro macro varibles to boundary nodes
for l=1:N
    NP=NODE{l};
    if NP{2}==0 % Interior nodes
        Rho_bc=inf;
        U_bc=[inf;inf];
        bc=Rho_bc;
        L=length(U_bc);
        bc(end+1:end+L,1)=U_bc;
        NP{21}=bc;
    elseif NP{2}<0 % Immersed boundary nodes
        if abs(NP{2})==4 || abs(NP{2})==73 % Stationary wall
            Rho_bc=inf;
            U_bc=[0;0];
            bc=Rho_bc;
            L=length(U_bc);
            bc(end+1:end+L,1)=U_bc;
            NP{21}=bc;
        elseif abs(NP{2})==5 || abs(NP{2})==73 % Moving wall
            error('The moving wall immersed boundary is not available!');
        else
            error('Some immersed boundary nodes are absent from macro bc variable definition!');
        end
    else % Outer boundary nodes
        %%%% Periodic and Fully developed, as interior nodes, nothing is known
        if NP{2}==1 || NP{2}==6
            Rho_bc=inf;
            U_bc=[inf;inf];
            bc=Rho_bc;
            L=length(U_bc);
            bc(end+1:end+L,1)=U_bc;
            NP{21}=bc;
            NS=NP{18};
            NUM=NODE{NS(1)};
            NDM=NODE{NS(2)};
            % Upstream node
            if (NUM{2}==1) || (NUM{2}==70) || (NUM{2}==71) || (NUM{2}==72)
                ;
            elseif NUM{2}==6
                ;
            elseif NUM{2}==76 || NUM{2}==77
                NUM{21}=bc;
                NODE{NUM{1}}=NUM;
            else
                error('Boundary node has incorrect numerical flag!');
            end
            % Downstream node
            if (NDM{2}==1) || (NDM{2}==70) || (NDM{2}==71) || (NDM{2}==72)
                ;
            elseif NDM{2}==6
                ;
            elseif NDM{2}==76 || NDM{2}==77
                NDM{21}=bc;
                NODE{NDM{1}}=NDM;
            else
                error('Boundary node has incorrect numerical flag!');
            end
            %%%% Velocity is known and constant, density is unknown
        elseif NP{2}==4
            Rho_bc=inf;
            U_bc=[0;0];
            bc=Rho_bc;
            L=length(U_bc);
            bc(end+1:end+L,1)=U_bc;
            NP{21}=bc;
            NS=NP{18};
            NUM=NODE{NS(1)};
            NDM=NODE{NS(2)};
            % Upstream node
            if NUM{2}==4
                ;
            elseif NUM{2}==71
                NUM{21}=bc;
                NODE{NUM{1}}=NUM;
            elseif NUM{2}==73
                if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                    Hydro_bc_corner(1:Lc,NUM{1})=bc;
                else
                    Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                end
            elseif NUM{2}==75
                Hydro_temp=NUM{21};
                Hydro_temp(2:Lc,1)=U_bc;
                NUM{21}=Hydro_temp;
                NODE{NUM{1}}=NUM;
            elseif NUM{2}==76
                ;
            else
                error('Boundary node has incorrect numerical flag!');
            end
            % Downstream node
            if NDM{2}==4
                ;
            elseif NDM{2}==71
                NDM{21}=bc;
                NODE{NDM{1}}=NDM;
            elseif NDM{2}==73
                if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                    Hydro_bc_corner(1:Lc,NDM{1})=bc;
                else
                    Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                end
            elseif NDM{2}==75
                Hydro_temp=NDM{21};
                Hydro_temp(2:Lc,1)=U_bc;
                NDM{21}=Hydro_temp;
                NODE{NDM{1}}=NDM;
            elseif NDM{2}==76
                ;
            else
                error('Boundary node has incorrect numerical flag!');
            end
        else
            if l<=N_I+N_L-2; %%%% Top
                %%%% Velocity is known but may varies at different locations, density is unknown
                if NP{2}==5
                    Rho_bc=inf;
                    U_bc=U_m(:,1); % First column for top surface
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==5
                        ;
                    elseif NUM{2}==71
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        if norm(Hydro_bc_corner(:,NUM{1}))==0 % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==5
                        ;
                    elseif NDM{2}==71
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==20
                    F=F_u_in(:,1); % First column for top surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_in(:,1); % First column for top surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,1);  % First column for top surface
                        U_max=U_in(:,1); % First column for top surface
                        U_x=U_ref(1,1)+sym_para(X1,X2,U_max(1,1),NC(1,1));
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_in(:,1); % First column for top surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,1);  % First column for top surface 
                        U_max=U_in(:,1); % First column for top surface
                        U_y=U_ref(2,1)+sym_para(X1,X2,U_max(2,1),NC(1,1));
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==20
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==20
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==30
                    F=F_u_out(:,1); % First column for top surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_out(:,1); % First column for top surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,1);  % First column for top surface
                        U_max=U_out(:,1); % First column for top surface
                        U_x=U_ref(1,1)+sym_para(X1,X2,U_max(1,1),NC(1,1));
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_out(:,1); % First column for top surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,1);  % First column for top surface
                        U_max=U_out(:,1); % First column for top surface
                        U_y=U_ref(2,1)+sym_para(X1,X2,U_max(2,1),NC(1,1));
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==30
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==30
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,1); % First column for top surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,1);  % First column for top surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                %%%% Density is known, velocity is unknown
                elseif NP{2}==21
                    F=F_rho_in(:,1); % First column for top surface
                    if F(:,1)==0
                        Rho_temp=Rho_in(:,1); % First column for top surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_in(:,1); % First column for top surface
                        Rho_temp=Rho_r+sym_para(X1,X2,Rho_max(:,1),NC(1,1)); % First column for top surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==21
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==21
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==31
                    F=F_rho_out(:,1); % First column for top surface
                    if F(:,1)==0
                        Rho_temp=Rho_out(:,1); % First column for top surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_out(:,1); % First column for top surface
                        Rho_temp=Rho_r+sym_para(X1,X2,Rho_max(:,1),NC(1,1)); % First column for top surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==31
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==31
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,1); % First column for top surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % First column for top surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                else
                    error('The current boundary node has incorrect numerical flag!');
                end
            elseif l>N_I+N_L-1 && l<=N_I+N_L-1+N_H-2; %%%% Right
                %%%% Velocity is known but may varies at different locations, density is unknown
                if NP{2}==5
                    Rho_bc=inf;
                    U_bc=U_m(:,2); % Second column for right surface
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==5
                        ;
                    elseif NUM{2}==71
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==5
                        ;
                    elseif NDM{2}==71
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==20
                    F=F_u_in(:,2); % Second column for right surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_in(:,2); % Second column for right surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,2);  % Second column for right surface
                        U_max=U_in(:,2); % Second column for right surface
                        U_x=U_ref(1,1)+sym_para(Y1,Y2,U_max(1,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_in(:,2); % Second column for right surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,2);  % First column for top surface
                        U_max=U_in(:,2); % First column for top surface
                        U_y=U_ref(2,1)+sym_para(Y1,Y2,U_max(2,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==20
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==20
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==30
                    F=F_u_out(:,2); % Second column for right surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_out(:,2); % Second column for right surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,2);  % Second column for right surface
                        U_max=U_out(:,2); % Second column for right surface
                        U_x=U_ref(1,1)+sym_para(Y1,Y2,U_max(1,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_out(:,2); % Second column for right surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,2);  % Second column for right surface
                        U_max=U_out(:,2); % Second column for top surface
                        U_y=U_ref(2,1)+sym_para(Y1,Y2,U_max(2,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==30
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==30
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,2); % Second column for right surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                %%%% Density is known, velocity is unknown
                elseif NP{2}==21
                    F=F_rho_in(:,2); % Second column for right surface
                    if F(:,1)==0
                        Rho_temp=Rho_in(:,2); % Second column for right surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_in(:,2); % Second column for right surface
                        Rho_temp=Rho_r+sym_para(Y1,Y2,Rho_max(:,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==21
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==21
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==31
                    F=F_rho_out(:,2); % Second column for right surface
                    if F(:,1)==0
                        Rho_temp=Rho_out(:,2); % Second column for right surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_out(:,2); % Second column for right surface
                        Rho_temp=Rho_r+sym_para(Y1,Y2,Rho_max(:,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==31
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==31
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,2); % Second column for right surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Second column for right surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                else
                    error('The current boundary node has incorrect numerical flag!');
                end
            elseif l>N_I+N_L-1+N_H-1 && l<=N_I+N_L-1+N_H-2+N_L-1; %%%% Bottom
                %%%% Velocity is known but may varies at different locations, density is unknown
                if NP{2}==5
                    Rho_bc=inf;
                    U_bc=U_m(:,3); % Third column for bottom surface
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==5
                        ;
                    elseif NUM{2}==71
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==5
                        ;
                    elseif NDM{2}==71
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==20
                    F=F_u_in(:,3); % Third column for bottom surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_in(:,3); % Third column for bottom surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,3);  % Third column for bottom surface
                        U_max=U_in(:,3); % Third column for bottom surface
                        U_x=U_ref(1,1)+sym_para(X1,X2,U_max(1,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_in(:,3); % Third column for bottom surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,3);  % Third column for bottom surface
                        U_max=U_in(:,3); % Third column for bottom surface
                        U_y=U_ref(2,1)+sym_para(X1,X2,U_max(2,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==20
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==20
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==30
                    F=F_u_out(:,3); % Third column for bottom surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_out(:,3); % Third column for bottom surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,3);  % Third column for bottom surface
                        U_max=U_out(:,3); % Third column for bottom surface
                        U_x=U_ref(1,1)+sym_para(X1,X2,U_max(1,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_out(:,3); % Third column for bottom surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,3);  % Third column for bottom surface
                        U_max=U_out(:,3); % Third column for bottom surface
                        U_y=U_ref(2,1)+sym_para(X1,X2,U_max(2,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==30
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==30
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,2);  % Second column for right surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,3); % Third column for bottom surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,3);  % Third column for bottom surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                %%%% Density is known, velocity is unknown
                elseif NP{2}==21
                    F=F_rho_in(:,3); % Third column for bottom surface
                    if F(:,1)==0
                        Rho_temp=Rho_in(:,3); % Third column for bottom surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_in(:,3); % Third column for bottom surface
                        Rho_temp=Rho_r+sym_para(X1,X2,Rho_max(:,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==21
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==21
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==31
                    F=F_rho_out(:,3); % Third column for bottom surface
                    if F(:,1)==0
                        Rho_temp=Rho_out(:,3); % Third column for bottom surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_out(:,3); % Third column for bottom surface
                        Rho_temp=Rho_r+sym_para(X1,X2,Rho_max(:,1),NC(1,1)); % Horizontal surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==31
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==31
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,3); % Third column for bottom surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Third column for bottom surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                else
                    error('The current boundary node has incorrect numerical flag!');
                end
            elseif  l>N_I+N_L-1+N_H-2+N_L && l<N  %%%% Left
                %%%% Velocity is known but may varies at different locations, density is unknown
                if NP{2}==5
                    Rho_bc=inf;
                    U_bc=U_m(:,4); % Fourth column for left surface
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==5
                        ;
                    elseif NUM{2}==71
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==5
                        ;
                    elseif NDM{2}==71
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        Hydro_temp(2:Lc,1)=U_bc;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==20
                    F=F_u_in(:,4); % Fourth column for left surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_in(:,4); % Fourth column for left surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,4);  % Fourth column for left surface
                        U_max=U_in(:,4); % Fourth column for left surface
                        U_x=U_ref(1,1)+sym_para(Y1,Y2,U_max(1,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_in(:,4); % Fourth column for left surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,4);  % Fourth column for left surface
                        U_max=U_in(:,4); % Fourth column for left surface
                        U_y=U_ref(2,1)+sym_para(Y1,Y2,U_max(2,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==20
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==20
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_in(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==30
                    F=F_u_out(:,4); % Fourth column for left surface
                    % X direction
                    if F(1,1)==0
                        U_temp=U_out(:,4); % Fourth column for left surface
                        U_x=U_temp(1,1);
                    elseif F(1,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(1,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,4);  % Fourth column for left surface
                        U_max=U_out(:,4); % Fourth column for left surface
                        U_x=U_ref(1,1)+sym_para(Y1,Y2,U_max(1,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    % Y direction
                    if F(2,1)==0
                        U_temp=U_out(:,4); % Fourth column for left surface
                        U_y=U_temp(2,1);
                    elseif F(2,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(2,1)==2
                        NC=NP{3};
                        U_ref=U_r(:,4);  % Fourth column for left surface
                        U_max=U_out(:,4); % Fourth column for left surface
                        U_y=U_ref(2,1)+sym_para(Y1,Y2,U_max(2,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=inf;
                    U_bc=[U_x;U_y];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==30
                        ;
                    elseif NUM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==30
                        ;
                    elseif NDM{2}==71
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==73
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(2:end,1)=[U_x;U_y];
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        % X direction
                        if F(1,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_x=U_temp(1,1);
                        elseif F(1,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(1,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_x=U_ref(1,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        % Y direction
                        if F(2,1)==0
                            U_temp=U_out(:,4); % Fourth column for left surface
                            U_y=U_temp(2,1);
                        elseif F(2,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(2,1)==2
                            U_ref=U_r(:,4);  % Fourth column for left surface
                            U_y=U_ref(2,1);
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(2:Lc,1)=[U_x;U_y];
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                %%%% Density is known, velocity is unknown
                elseif NP{2}==21
                    F=F_rho_in(:,4); % Fourth column for left surface
                    if F(:,1)==0
                        Rho_temp=Rho_in(:,4); % Fourth column for left surface
                    elseif F(:,1)==1
                        NC=NP{3};
                        Rho1=Rho_in(:,3);
                        Rho2=Rho_in(:,1);
                        a=(Rho1-Rho2)/(Y1-Y2);
                        b=Rho1-Y1*a;
                        Rho_temp=a*NC(2,1)+b;
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_in(:,4); % Fourth column for left surface
                        Rho_temp=Rho_r+sym_para(Y1,Y2,Rho_max(:,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==21
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            Rho_temp=Rho_in(:,3);
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==21
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_in(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            Rho_temp=Rho_in(:,1);
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                elseif NP{2}==31
                    F=F_rho_out(:,4); % Fourth column for left surface
                    if F(:,1)==0
                        Rho_temp=Rho_out(:,4); % Fourth column for left surface
                    elseif F(:,1)==1
                        error('Temporarilly unavailable!');
                    elseif F(:,1)==2
                        NC=NP{3};
                        Rho_max=Rho_out(:,4); % Fourth column for left surface
                        Rho_temp=Rho_r+sym_para(Y1,Y2,Rho_max(:,1),NC(2,1)); % Vertical surface
                    else
                        error('The flag for velocity profile is invalid or unavailable!');
                    end
                    Rho_bc=Rho_temp;
                    U_bc=[Inf;Inf];
                    bc=Rho_bc;
                    L=length(U_bc);
                    bc(end+1:end+L,1)=U_bc;
                    NP{21}=bc;
                    NS=NP{18};
                    NUM=NODE{NS(1)};
                    NDM=NODE{NS(2)};
                    % Upstream node
                    if NUM{2}==31
                        ;
                    elseif NUM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NUM{21}=bc;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NUM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NUM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NUM{1})=bc;
                        end
                    elseif NUM{2}==75
                        Hydro_temp=NUM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NUM{21}=Hydro_temp;
                        NODE{NUM{1}}=NUM;
                    elseif NUM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                    % Downstream node
                    if NDM{2}==31
                        ;
                    elseif NDM{2}==72
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r; % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        NDM{21}=bc;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==74
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        bc(1,1)=Rho_temp;
                        if norm(Hydro_bc_corner(:,NDM{1}))==0; % Still empty
                            Hydro_bc_corner(1:Lc,NDM{1})=bc;
                        else
                            Hydro_bc_corner(Lc+1:2*Lc,NDM{1})=bc;
                        end
                    elseif NDM{2}==75
                        Hydro_temp=NDM{21};
                        if F(:,1)==0
                            Rho_temp=Rho_out(:,4); % Fourth column for left surface
                        elseif F(:,1)==1
                            error('Temporarilly unavailable!');
                        elseif F(:,1)==2
                            Rho_temp=Rho_r;  % Fourth column for left surface
                        else
                            error('The flag for velocity profile is invalid or unavailable!');
                        end
                        Hydro_temp(1,1)=Rho_temp;
                        NDM{21}=Hydro_temp;
                        NODE{NDM{1}}=NDM;
                    elseif NDM{2}==76
                        ;
                    else
                        error('Boundary node has incorrect numerical flag!');
                    end
                else
                    error('The current boundary node has incorrect numerical flag!');
                end
            elseif l==N_I+N_L-1 % Top right
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
            elseif l==N_I+N_L-1+N_H-1 % Bottom right
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
            elseif l==N_I+N_L-1+N_H-2+N_L % Bottom left
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
            elseif l==N % Top left
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
            else
                error('Some outer boundary nodes are missing!');
            end
        end
    end
    NODE{l}=NP;
end
% Fill the data into corner nodes if their flags are 73 or 74
for l=1:N
    if norm(Hydro_bc_corner(:,l))~=0
        Hydro1=Hydro_bc_corner(1:Lc,l);
        Hydro2=Hydro_bc_corner(Lc+1:2*Lc,l);
        rho1=Hydro1(1,1);
        rho2=Hydro2(1,1);
        u1=Hydro1(2,1);
        u2=Hydro2(2,1);
        v1=Hydro1(3,1);
        v2=Hydro2(3,1);
        ND=NODE{l};
        if l==N_I+N_L-1 % Top right
            if FCO(1)==0
                rho=min(rho1,rho2);
                u=min(u1,u2);
                v=min(v1,v2);
            elseif FCO(1)==1
                rho=max(rho1,rho2);
                u=max(u1,u2);
                v=max(v1,v2);
            elseif FCO(1)==2
                rho=0.5*(rho1+rho2);
                u=0.5*(u1+u2);
                v=0.5*(v1+v2);
            else
                error('Wrong flag!');
            end
        elseif l==N_I+N_L-1+N_H-1 % Bottom right
            if FCO(2)==0
                rho=min(rho1,rho2);
                u=min(u1,u2);
                v=min(v1,v2);
            elseif FCO(2)==1
                rho=max(rho1,rho2);
                u=max(u1,u2);
                v=max(v1,v2);
            elseif FCO(2)==2
                rho=0.5*(rho1+rho2);
                u=0.5*(u1+u2);
                v=0.5*(v1+v2);
            else
                error('Wrong flag!');
            end
        elseif l==N_I+N_L-1+N_H-2+N_L % Bottom left
            if FCO(3)==0
                rho=min(rho1,rho2);
                u=min(u1,u2);
                v=min(v1,v2);
            elseif FCO(3)==1
                rho=max(rho1,rho2);
                u=max(u1,u2);
                v=max(v1,v2);
            elseif FCO(3)==2
                rho=0.5*(rho1+rho2);
                u=0.5*(u1+u2);
                v=0.5*(v1+v2);
            else
                error('Wrong flag!');
            end
        elseif l==N % Top left
            if FCO(4)==0
                rho=min(rho1,rho2);
                u=min(u1,u2);
                v=min(v1,v2);
            elseif FCO(4)==1
                rho=max(rho1,rho2);
                u=max(u1,u2);
                v=max(v1,v2);
            elseif FCO(4)==2
                rho=0.5*(rho1+rho2);
                u=0.5*(u1+u2);
                v=0.5*(v1+v2);
            else
                error('Wrong flag!');
            end
        else
            error('The node is not corner node!');
        end
        Hydro=[rho;u;v];
        ND{21}=Hydro;
        NODE{l}=ND;
    end
end


%% Thermal
T_top=T_bc(1,1);
T_right=T_bc(1,2);
T_bottom=T_bc(1,3);
T_left=T_bc(1,4);
T_imm=T_bc(1,5);
% Check positivity of temperature
if ((T_top<0 || T_bottom<0) || (T_left<0 || T_right<0)) || T_imm<0
    error('Temperature cannot be negative!');
end
% Temporary container for corner nodes
Therm_bc_corner=zeros(2,N);
% Write thermal macro varibles to boundary nodes
for l=1:N
    NP=NODE{l};
    if NP{2}==0 % Interior nodes
        bc=NP{21};
        bc(4)=inf;
        NP{21}=bc;
    elseif NP{2}<0 % Immersed boundary nodes
        bc=NP{21};
        bc(4)=T_imm;
        NP{21}=bc;
    else % Outer boundary nodes
        %%%% Periodic and Fully developed, as interior nodes, nothing is known
        if NP{2}==1 || NP{2}==6
            bc=NP{21};
            bc(4)=inf;
            NP{21}=bc;
            NS=NP{18};
            NUM=NODE{NS(1)};
            NDM=NODE{NS(2)};
            % Upstream node
            if (NUM{2}==1) || (NUM{2}==70) || (NUM{2}==71) || (NUM{2}==72)
                ;
            elseif NUM{2}==6
                ;
            elseif NUM{2}==76 || NUM{2}==77
                bc=NUM{21};
                bc(4)=inf;
                NUM{21}=bc;
                NODE{NUM{1}}=NUM;
            else
                error('Boundary node has incorrect numerical flag!');
            end
            % Downstream node
            if (NDM{2}==1) || (NDM{2}==70) || (NDM{2}==71) || (NDM{2}==72)
                ;
            elseif NDM{2}==6
                ;
            elseif NDM{2}==76 || NDM{2}==77
                bc=NDM{21};
                bc(4)=inf;
                NDM{21}=bc;
                NODE{NDM{1}}=NDM;
            else
                error('Boundary node has incorrect numerical flag!');
            end
        else
            if l<=N_I+N_L-2; %%%% Top
                bc=NP{21};
                bc(4)=T_top;
                NP{21}=bc;
                NS=NP{18};
                %% Up and down stream nodes
                NUM=NODE{NS(1)};
                NDM=NODE{NS(2)};
                % Upstream node
                if NUM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NUM{2}==71 || NUM{2}==72
                    bc=NUM{21};
                    bc(4)=T_top;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==73 || (NUM{2}==74 || NUM{2}==75)
                    if norm(Therm_bc_corner(:,NUM{1}))==0; % Still empty
                        Therm_bc_corner(1,NUM{1})=T_top;
                    else
                        Therm_bc_corner(2,NUM{1})=T_top;
                    end
                elseif NUM{2}==76
                    bc=NUM{21};
                    bc(4)=inf;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NUM{2}==1 || NUM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
                % Downstream node
                if NDM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NDM{2}==71 || NDM{2}==72
                    bc=NDM{21};
                    bc(4)=T_top;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==73 || (NDM{2}==74 || NDM{2}==75)
                    if norm(Therm_bc_corner(:,NDM{1}))==0; % Still empty
                        Therm_bc_corner(1,NDM{1})=T_top;
                    else
                        Therm_bc_corner(2,NDM{1})=T_top;
                    end
                elseif NDM{2}==76
                    bc=NDM{21};
                    bc(4)=inf;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NDM{2}==1 || NDM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
            elseif l>N_I+N_L-1 && l<=N_I+N_L-1+N_H-2; %%%% Right
                bc=NP{21};
                bc(4)=T_right;
                NP{21}=bc;
                NS=NP{18};
                %% Up and down stream nodes
                NUM=NODE{NS(1)};
                NDM=NODE{NS(2)};
                % Upstream node
                if NUM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NUM{2}==71 || NUM{2}==72
                    bc=NUM{21};
                    bc(4)=T_right;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==73 || (NUM{2}==74 || NUM{2}==75)
                    if norm(Therm_bc_corner(:,NUM{1}))==0; % Still empty
                        Therm_bc_corner(1,NUM{1})=T_right;
                    else
                        Therm_bc_corner(2,NUM{1})=T_right;
                    end
                elseif NUM{2}==76
                    bc=NUM{21};
                    bc(4)=inf;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NUM{2}==1 || NUM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
                % Downstream node
                if NDM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NDM{2}==71 || NDM{2}==72
                    bc=NDM{21};
                    bc(4)=T_right;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==73 || (NDM{2}==74 || NDM{2}==75)
                    if norm(Therm_bc_corner(:,NDM{1}))==0; % Still empty
                        Therm_bc_corner(1,NDM{1})=T_right;
                    else
                        Therm_bc_corner(2,NDM{1})=T_right;
                    end
                elseif NDM{2}==76
                    bc=NDM{21};
                    bc(4)=inf;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NDM{2}==1 || NDM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
            elseif l>N_I+N_L-1+N_H-1 && l<=N_I+N_L-1+N_H-2+N_L-1; %%%% Bottom
                bc=NP{21};
                bc(4)=T_bottom;
                NP{21}=bc;
                NS=NP{18};
                %% Up and down stream nodes
                NUM=NODE{NS(1)};
                NDM=NODE{NS(2)};
                % Upstream node
                if NUM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NUM{2}==71 || NUM{2}==72
                    bc=NUM{21};
                    bc(4)=T_bottom;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==73 || (NUM{2}==74 || NUM{2}==75)
                    if norm(Therm_bc_corner(:,NUM{1}))==0; % Still empty
                        Therm_bc_corner(1,NUM{1})=T_bottom;
                    else
                        Therm_bc_corner(2,NUM{1})=T_bottom;
                    end
                elseif NUM{2}==76
                    bc=NUM{21};
                    bc(4)=inf;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NUM{2}==1 || NUM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
                % Downstream node
                if NDM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NDM{2}==71 || NDM{2}==72
                    bc=NDM{21};
                    bc(4)=T_bottom;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==73 || (NDM{2}==74 || NDM{2}==75)
                    if norm(Therm_bc_corner(:,NDM{1}))==0; % Still empty
                        Therm_bc_corner(1,NDM{1})=T_bottom;
                    else
                        Therm_bc_corner(2,NDM{1})=T_bottom;
                    end
                elseif NDM{2}==76
                    bc=NDM{21};
                    bc(4)=inf;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NDM{2}==1 || NDM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
            elseif  l>N_I+N_L-1+N_H-2+N_L && l<N  %%%% Left
                bc=NP{21};
                bc(4)=T_left;
                NP{21}=bc;
                NS=NP{18};
                %% Up and down stream nodes
                NUM=NODE{NS(1)};
                NDM=NODE{NS(2)};
                % Upstream node
                if NUM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NUM{2}==71 || NUM{2}==72
                    bc=NUM{21};
                    bc(4)=T_left;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==73 || (NUM{2}==74 || NUM{2}==75)
                    if norm(Therm_bc_corner(:,NUM{1}))==0; % Still empty
                        Therm_bc_corner(1,NUM{1})=T_left;
                    else
                        Therm_bc_corner(2,NUM{1})=T_left;
                    end
                elseif NUM{2}==76
                    bc=NUM{21};
                    bc(4)=inf;
                    NUM{21}=bc;
                    NODE{NUM{1}}=NUM;
                elseif NUM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NUM{2}==1 || NUM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
                % Downstream node
                if NDM{2}==70
                    error('Boundary node has incorrect numerical flag!');
                elseif NDM{2}==71 || NDM{2}==72
                    bc=NDM{21};
                    bc(4)=T_left;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==73 || (NDM{2}==74 || NDM{2}==75)
                    if norm(Therm_bc_corner(:,NDM{1}))==0; % Still empty
                        Therm_bc_corner(1,NDM{1})=T_left;
                    else
                        Therm_bc_corner(2,NDM{1})=T_left;
                    end
                elseif NDM{2}==76
                    bc=NDM{21};
                    bc(4)=inf;
                    NDM{21}=bc;
                    NODE{NDM{1}}=NDM;
                elseif NDM{2}==77
                    error('Boundary node has incorrect numerical flag!');
                else
                    if NDM{2}==1 || NDM{2}==6
                        error('Boundary node has incorrect numerical flag!');
                    else
                        ;
                    end
                end
            elseif l==N_I+N_L-1 % Top right
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
                if NP{2}==71
                    ;
                else
                    bc=NP{21};
                    bc(4)=T_right;
                    NP{21}=bc;
                end
            elseif l==N_I+N_L-1+N_H-1 % Bottom right
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
                if NP{2}==71
                    ;
                else
                    bc=NP{21};
                    bc(4)=T_right;
                    NP{21}=bc;
                end
            elseif l==N_I+N_L-1+N_H-2+N_L % Bottom left
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
                if NP{2}==71
                    ;
                else
                    bc=NP{21};
                    bc(4)=T_left;
                    NP{21}=bc;
                end
            elseif l==N % Top left
                if NP{2}<70
                    error('The current corner node has incorrect numerical flag or current node is not at corners!');
                end
                if NP{2}==71
                    ;
                else
                    bc=NP{21};
                    bc(4)=T_left;
                    NP{21}=bc;
                end
            else
                error('Some outer boundary nodes are missing!');
            end
        end
    end
    NODE{l}=NP;
end
% Fill the data into corner nodes if their flags are 73, 74 or 75
for l=1:N
    if norm(Therm_bc_corner(:,l))~=0
        T1=Therm_bc_corner(1,l);
        T2=Therm_bc_corner(2,l);
        ND=NODE{l};
        if l==N_I+N_L-1 % Top right
            if FCO(1)==0
                T=min(T1,T2);
            elseif FCO(1)==1
                T=max(T1,T2);
            elseif FCO(1)==2
                T=0.5*(T1+T2);
            else
                error('Wrong flag!');
            end
        elseif l==N_I+N_L-1+N_H-1 % Bottom right
            if FCO(2)==0
                T=min(T1,T2);
            elseif FCO(2)==1
                T=max(T1,T2);
            elseif FCO(2)==2
                T=0.5*(T1+T2);
            else
                error('Wrong flag!');
            end
        elseif l==N_I+N_L-1+N_H-2+N_L % Bottom left
            if FCO(3)==0
                T=min(T1,T2);
            elseif FCO(3)==1
                T=max(T1,T2);
            elseif FCO(3)==2
                T=0.5*(T1+T2);
            else
                error('Wrong flag!');
            end
        elseif l==N % Top left
            if FCO(4)==0
                T=min(T1,T2);
            elseif FCO(4)==1
                T=max(T1,T2);
            elseif FCO(4)==2
                T=0.5*(T1+T2);
            else
                error('Wrong flag!');
            end
        else
            error('The node is not corner node!');
        end
        bc=ND{21};
        bc(4)=T;
        ND{21}=bc;
        NODE{l}=ND;
    end
end