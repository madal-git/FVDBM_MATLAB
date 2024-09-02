function [CELL,NODE,FACE,FPDC]=bc(CELL,M,NODE,N,FACE,O,N_I,N_I_N,N_L,N_H,X1,X2,Y1,Y2,h,dx,dy,FM,V1,V2,top,right,bottom,left,immersed)
% function [CELL,NODE,FACE]=bc(CELL,M,NODE,N,FACE,O,N_I,N_I_N,N_L,N_H,X1,X2,Y1,Y2,h,dx,dy,FM,V1,V2,top,right,bottom,left,immersed)
% defines the physical boundary conditions of the model and writes those information
% into the triangle data structure CELL and nodal data structure NODE.
% CELL is cell structure data contains all cells and their information
% M is the total number of triangles
% NODE is cell structure data contains all nodes and their information
% N is the total number of nodes
% FACE is cell structure data contains all faces and their information
% O is the total number of nodes
% N_I is the total number of interior nodes (including the immersed boundary nodes)
% N_I_N is the total number of immersed boundary nodes
% N_L is the number of nodes on the length face
% N_H is the number of nodes on the height face
% the boundary nodes
% X1, X2, Y1 and Y2 are bounds of mesh domain
% h is the node spacing on outer boundaries for IRT mesh
% dx, dy is the node spacing on outer boundaries for random mesh.
% FM  is the mesh type flag, FM=0---IRT mesh;FM=1---random mesh
% FPDC is the flag for periodic boundary conditions. FPDC=0---No periodic
% boundaries; FPDC=1---Only left & right boundaries are periodic;
% FPDC=2---Only top & bottom boundaries are periodic; FPDC=3---All
% boundaries are periodic
%%%%%%% The following variables are string and all case-sensitive %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% top is string type of boundary condition on the top surface of the rectangular domain
% right is string type of boundary condition on the right surface of the rectangular domain
% bottom is string type of boundary condition on the bottom surface of the rectangular domain
% left is string type of boundary condition on the left surface of the rectangular domain
% immersed is string type of boundary condition on the immersed surface inside the rectangular domain

%%%% There are 6 types of boundary conditions as of now. 'Periodic',
%%%% 'Inlet', 'Outlet', 'Stationary Wall', 'Moving Wall' and 'Fully Developed'. 

%%%% Definitions and corresponding numeric flag of different boundary condition %%%% 
% For the nodes
% ND=NODE{l}
% ND{2}  Store the numeric flag of boundary condition of current node
% ND{2}=0         Interior
% ND{2}=1         Periodic  (has to be paired)
% ND{2}=20        Velocity Inlet
% ND{2}=21        Density Inlet
% ND{2}=30        Velocity Outlet
% ND{2}=31        Density Outlet
% abs(ND{2})=4    Stationary Wall
% abs(ND{2})=5    Moving Wall
% ND{2}=6         Fully Developed

% For the cells
% C=CELL{l}
% C{10},C{11} and C{12}  Store the numeric flag of boundary condition of three nodes of current triangle
% C{19},C{20} and C{21}  Store the numeric flag of boundary condition of three faces of current triangle
% C{10, 11 or 12}=0          Interior
% C{10, 11 or 12}=1          Periodic  (has to be paired)
% C{10, 11 or 12}=20         Velocity Inlet
% C{10, 11 or 12}=21         Density Inlet
% C{10, 11 or 12}=30         Velocity Outlet
% C{10, 11 or 12}=31         Density Outlet
% abs(C{10, 11 or 12})=4     Stationary Wall
% abs(C{10, 11 or 12})=5     Moving Wall
% C{10, 11 or 12}=6          Fully Developed

% C{19, 20 or 21}=0          Interior
% C{19, 20 or 21}=1          Periodic  (has to be paired)
% C{19, 20 or 21}=20         Velocity Inlet
% C{19, 20 or 21}=21         Density Inlet
% C{19, 20 or 21}=30         Velocity Outlet
% C{19, 20 or 21}=31         Density Outlet
% abs(C{19, 20 or 21})=4     Stationary Wall
% abs(C{19, 20 or 21})=5     Moving Wall
% C{19, 20 or 21}=6          Fully Developed

% For the faces
% FC=FACE{l}
% FC{2}  Store the numeric flag of boundary condition of current node
% FC{2}=0         Interior
% FC{2}=1         Periodic  (has to be paired)
% FC{2}=20        Velocity Inlet
% FC{2}=21        Density Inlet
% FC{2}=30        Velocity Outlet
% FC{2}=31        Density Outlet
% abs(FC{2})=4    Stationary Wall
% abs(FC{2})=5    Moving Wall
% FC{2}=6         Fully Developed

q1=length(V1(1,:));
q2=length(V2(1,:));

P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size

TEMP=10; % Has to be integer but NOT one
%% Check pair of 'Periodic' boundary condition %%%%
pp=0; % Number of periodic boundary condition
if strcmp(top,'Periodic')==1
    pp=pp+1;
end
if strcmp(right,'Periodic')==1
    pp=pp+1;
end
if strcmp(bottom,'Periodic')==1
    pp=pp+1;
end
if strcmp(left,'Periodic')==1
    pp=pp+1;
end
% Check pair and face-to-face of 'Periodic' boundary condition %%%%
if mod(pp,2)~=0
    error('The Periodic boundary conditions are not paired!');
else
    if (strcmp(top,'Periodic')==1 && strcmp(bottom,'Periodic')~=1) || (strcmp(top,'Periodic')~=1 && strcmp(bottom,'Periodic')==1)
        error('The exsisting periodic pair is not face-to-face!');
    elseif (strcmp(left,'Periodic')==1 && strcmp(right,'Periodic')~=1) || (strcmp(left,'Periodic')~=1 && strcmp(right,'Periodic')==1)
        error('The exsisting periodic pair is not face-to-face!');
    else
        ;
    end
end
%% Output the global flag for periodic boundary conditions
if pp==0
    FPDC=0;
elseif pp==2
    if (strcmp(top,'Periodic')~=1 && strcmp(bottom,'Periodic')~=1) && (strcmp(left,'Periodic')==1 && strcmp(right,'Periodic')==1)
        FPDC=1;
    elseif (strcmp(top,'Periodic')==1 && strcmp(bottom,'Periodic')==1) && (strcmp(left,'Periodic')~=1 && strcmp(right,'Periodic')~=1)
        FPDC=2;
    else
        error('Logic error!');
    end
elseif pp==4
    FPDC=3;
else
    error('Logic error!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary flags are added saved to NODE data structure
%% Nodes on the plat boundaries
for l=1:N;
    NP=NODE{l};
    if NP{2}==0 %Interior nodes
        if NP{1}<=N_I_N || NP{1}>N_I
            error('The boundary identifier for interior node is not ZERO!');
        end
    elseif NP{2}>0 %Outer boundary nodes
        if NP{1}<=N_I
            error('The boundary identifier for outer boundary node is not POSITIVE!');
        end
        if l<=N_I+N_L-2; %%%% Top
            if strcmp(top,'Periodic')==1
                NP{2}=1;
            elseif strcmp(top,'Velocity Inlet')==1
                NP{2}=20;
            elseif strcmp(top,'Pressure Inlet')==1
                NP{2}=21;
            elseif strcmp(top,'Velocity Outlet')==1
                NP{2}=30;
            elseif strcmp(top,'Pressure Outlet')==1
                NP{2}=31;
            elseif strcmp(top,'Stationary Wall')==1
                NP{2}=4;
            elseif strcmp(top,'Moving Wall')==1
                NP{2}=5;
            elseif strcmp(top,'Fully Developed')==1
                NP{2}=6;
            else
                error('Please check spelling and input the correct boundary type!');
            end
        elseif l>N_I+N_L-1 && l<=N_I+N_L-1+N_H-2; %%%% Right
            if strcmp(right,'Periodic')==1
                NP{2}=1;
            elseif strcmp(right,'Velocity Inlet')==1
                NP{2}=20;
            elseif strcmp(right,'Pressure Inlet')==1
                NP{2}=21;
            elseif strcmp(right,'Velocity Outlet')==1
                NP{2}=30;
            elseif strcmp(right,'Pressure Outlet')==1
                NP{2}=31;
            elseif strcmp(right,'Stationary Wall')==1
                NP{2}=4;
            elseif strcmp(right,'Moving Wall')==1
                NP{2}=5;
            elseif strcmp(right,'Fully Developed')==1
                NP{2}=6;
            else
                error('Please check spelling and input the correct boundary type!');
            end
        elseif l>N_I+N_L-1+N_H-1 && l<=N_I+N_L-1+N_H-2+N_L-1; %%%% Bottom
            if strcmp(bottom,'Periodic')==1
                NP{2}=1;
            elseif strcmp(bottom,'Velocity Inlet')==1
                NP{2}=20;
            elseif strcmp(bottom,'Pressure Inlet')==1
                NP{2}=21;
            elseif strcmp(bottom,'Velocity Outlet')==1
                NP{2}=30;
            elseif strcmp(bottom,'Pressure Outlet')==1
                NP{2}=31;
            elseif strcmp(bottom,'Stationary Wall')==1
                NP{2}=4;
            elseif strcmp(bottom,'Moving Wall')==1
                NP{2}=5;
            elseif strcmp(bottom,'Fully Developed')==1
                NP{2}=6;
            else
                error('Please check spelling and input the correct boundary type!');
            end
        elseif  l>N_I+N_L-1+N_H-2+N_L && l<N  %%%% Left
            if strcmp(left,'Periodic')==1
                NP{2}=1;
            elseif strcmp(left,'Velocity Inlet')==1
                NP{2}=20;
            elseif strcmp(left,'Pressure Inlet')==1
                NP{2}=21;
            elseif strcmp(left,'Velocity Outlet')==1
                NP{2}=30;
            elseif strcmp(left,'Pressure Outlet')==1
                NP{2}=31;
            elseif strcmp(left,'Stationary Wall')==1
                NP{2}=4;
            elseif strcmp(left,'Moving Wall')==1
                NP{2}=5;
            elseif strcmp(left,'Fully Developed')==1
                NP{2}=6;
            else
                error('Please check spelling and input the correct boundary type!');
            end
        elseif l==N_I+N_L-1
            NP{2}=TEMP; % Temporary flag
        elseif l==N_I+N_L-1+N_H-1
            NP{2}=TEMP; % Temporary flag
        elseif l==N_I+N_L-1+N_H-2+N_L
            NP{2}=TEMP; % Temporary flag
        elseif l==N
            NP{2}=TEMP; % Temporary flag
        else
            error('Some nodes are absent from BC definition!');
        end
%     elseif NP{2}==1 %Outer boundary nodes
%         if l<=N_I+N_L-2; %%%% Top
%             if strcmp(top,'Periodic')==1
%                 NP{2}=NP{2}*1;
%             elseif strcmp(top,'Velocity Inlet')==1
%                 NP{2}=NP{2}*20;
%             elseif strcmp(top,'Pressure Inlet')==1
%                 NP{2}=NP{2}*21;
%             elseif strcmp(top,'Velocity Outlet')==1
%                 NP{2}=NP{2}*30;
%             elseif strcmp(top,'Pressure Outlet')==1
%                 NP{2}=NP{2}*31;
%             elseif strcmp(top,'Stationary Wall')==1
%                 NP{2}=NP{2}*4;
%             elseif strcmp(top,'Moving Wall')==1
%                 NP{2}=NP{2}*5;
%             elseif strcmp(top,'Fully Developed')==1
%                 NP{2}=NP{2}*6;
%             else
%                 error('Please check spelling and input the correct boundary type!');
%             end
%         elseif l>N_I+N_L-1 && l<=N_I+N_L-1+N_H-2; %%%% Right
%             if strcmp(right,'Periodic')==1
%                 NP{2}=NP{2}*1;
%             elseif strcmp(right,'Velocity Inlet')==1
%                 NP{2}=NP{2}*20;
%             elseif strcmp(right,'Pressure Inlet')==1
%                 NP{2}=NP{2}*21;
%             elseif strcmp(right,'Velocity Outlet')==1
%                 NP{2}=NP{2}*30;
%             elseif strcmp(right,'Pressure Outlet')==1
%                 NP{2}=NP{2}*31;
%             elseif strcmp(right,'Stationary Wall')==1
%                 NP{2}=NP{2}*4;
%             elseif strcmp(right,'Moving Wall')==1
%                 NP{2}=NP{2}*5;
%             elseif strcmp(right,'Fully Developed')==1
%                 NP{2}=NP{2}*6;
%             else
%                 error('Please check spelling and input the correct boundary type!');
%             end
%         elseif l>N_I+N_L-1+N_H-1 && l<=N_I+N_L-1+N_H-2+N_L-1; %%%% Bottom
%             if strcmp(bottom,'Periodic')==1
%                 NP{2}=NP{2}*1;
%             elseif strcmp(bottom,'Velocity Inlet')==1
%                 NP{2}=NP{2}*20;
%             elseif strcmp(bottom,'Pressure Inlet')==1
%                 NP{2}=NP{2}*21;
%             elseif strcmp(bottom,'Velocity Outlet')==1
%                 NP{2}=NP{2}*30;
%             elseif strcmp(bottom,'Pressure Outlet')==1
%                 NP{2}=NP{2}*31;
%             elseif strcmp(bottom,'Stationary Wall')==1
%                 NP{2}=NP{2}*4;
%             elseif strcmp(bottom,'Moving Wall')==1
%                 NP{2}=NP{2}*5;
%             elseif strcmp(bottom,'Fully Developed')==1
%                 NP{2}=NP{2}*6;
%             else
%                 error('Please check spelling and input the correct boundary type!');
%             end
%         elseif  l>N_I+N_L-1+N_H-2+N_L && l<N  %%%% Left
%             if strcmp(left,'Periodic')==1
%                 NP{2}=NP{2}*1;
%             elseif strcmp(left,'Velocity Inlet')==1
%                 NP{2}=NP{2}*20;
%             elseif strcmp(left,'Pressure Inlet')==1
%                 NP{2}=NP{2}*21;
%             elseif strcmp(left,'Velocity Outlet')==1
%                 NP{2}=NP{2}*30;
%             elseif strcmp(left,'Pressure Outlet')==1
%                 NP{2}=NP{2}*31;
%             elseif strcmp(left,'Stationary Wall')==1
%                 NP{2}=NP{2}*4;
%             elseif strcmp(left,'Moving Wall')==1
%                 NP{2}=NP{2}*5;
%             elseif strcmp(left,'Fully Developed')==1
%                 NP{2}=NP{2}*6;
%             else
%                 error('Please check spelling and input the correct boundary type!');
%             end
%         elseif l==N_I+N_L-1
%             NP{2}=NP{2}*TEMP; % Temporary flag
%         elseif l==N_I+N_L-1+N_H-1
%             NP{2}=NP{2}*TEMP; % Temporary flag
%         elseif l==N_I+N_L-1+N_H-2+N_L
%             NP{2}=NP{2}*TEMP; % Temporary flag
%         elseif l==N
%             NP{2}=NP{2}*TEMP; % Temporary flag
%         else
%             error('Some nodes are absent from BC definition!');
%         end
    elseif NP{2}<0 %Immersed boundary nodes
        if NP{1}>N_I_N
            error('The boundary identifier for immersed boundary node is not NEGATIVE!');
        end
        if strcmp(immersed,'Periodic')==1
            error('Immersed boundary cannot be periodic!');
        elseif strcmp(immersed,'Velocity Inlet')==1
            error('Temporarily not available!');
        elseif strcmp(immersed,'Pressure Inlet')==1
            error('Temporarily not available!');
        elseif strcmp(immersed,'Velocity Outlet')==1
            error('Temporarily not available!');
        elseif strcmp(immersed,'Pressure Outlet')==1
            error('Temporarily not available!');
        elseif strcmp(immersed,'Stationary Wall')==1
            NP{2}=-4;
        elseif strcmp(immersed,'Moving Wall')==1
            error('Temporarily not available!');
        elseif strcmp(immersed,'Fully Developed')==1
            error('Immersed boundary cannot be fully developed!');
        else
            error('Temporarily not available!');
        end
%     elseif NP{2}==-1 %Immersed boundary nodes
%         if strcmp(immersed,'Periodic')==1
%             error('Immersed boundary cannot be periodic!');
%         elseif strcmp(immersed,'Velocity Inlet')==1
%             error('Temporarily not available!');
%         elseif strcmp(immersed,'Pressure Inlet')==1
%             error('Temporarily not available!');
%         elseif strcmp(immersed,'Velocity Outlet')==1
%             error('Temporarily not available!');
%         elseif strcmp(immersed,'Pressure Outlet')==1
%             error('Temporarily not available!');
%         elseif strcmp(immersed,'Stationary Wall')==1
%             NP{2}=NP{2}*4;
%         elseif strcmp(immersed,'Moving Wall')==1
%             error('Temporarily not available!');
%         elseif strcmp(immersed,'Fully Developed')==1
%             error('Immersed boundary cannot be fully developed!');
%         else
%             error('Temporarily not available!');
%         end
    else
        error('The initial boundary flags are misslabled! Go back to meshing tools to check!');
    end
    NODE{l}=NP;
end
%%% Corner nodes
for l=1:N
    NP=NODE{l};
    if abs(NP{2})==TEMP % Corner nodes
        if l==N_I+N_L-1 || l==N_I+N_L-1+N_H-1 || l==N_I+N_L-1+N_H-2+N_L || l==N %%%% Top right corner, Bottom right  corner, Bottom left  corner, Top left corner. Can be replaced with if NP{2}>=0
            NS=NP{18};
            NU=NODE{NS(1)};
            ND=NODE{NS(2)};
            if abs(NU{2})==6 || abs(ND{2})==6 % At least one side of the corner node is fully developed
                if abs(NU{2})==6 && abs(ND{2})==6
                    NP{2}=(NP{2}/TEMP)*77; % Both sides are fully-developed
                else
                    NP{2}=(NP{2}/TEMP)*76; % Only one side is fully-developed
                end
            elseif abs(NU{2})==1 || abs(ND{2})==1
                if abs(NU{2})==1 && abs(ND{2})==1
                    NP{2}=(NP{2}/TEMP)*70; % Both sides are periodic
                elseif (abs(NU{2})==4 || abs(NU{2})==5 || abs(NU{2})==20 || abs(NU{2})==30) || (abs(ND{2})==4 || abs(ND{2})==5 || abs(ND{2})==20 || abs(ND{2})==30)
                    % One side is periodic, the other side is velocity
                    NP{2}=(NP{2}/TEMP)*71;
                elseif (abs(NU{2})==21 || abs(NU{2})==31) || (abs(ND{2})==21 || abs(ND{2})==31)
                    % One side is periodic, the other side is density
                    NP{2}=(NP{2}/TEMP)*72;
                else
                    error('Some corner nodes boardering periodic BC are absent from BC definition!')
                end
            elseif (abs(NU{2})==4 || abs(NU{2})==5 || abs(NU{2})==20 || abs(NU{2})==30) && (abs(ND{2})==4 || abs(ND{2})==5 || abs(ND{2})==20 || abs(ND{2})==30)
                % Both sides are velocity
                NP{2}=(NP{2}/TEMP)*73;
            elseif (abs(NU{2})==21 || abs(NU{2})==31) && (abs(ND{2})==21 || abs(ND{2})==31)
                % Both sides are density
                NP{2}=(NP{2}/TEMP)*74;
            else % One side is velocity, the other is density
                NP{2}=(NP{2}/TEMP)*75;
            end
        else
            error('Some corner nodes are absent from BC definition!');
        end
    end
    NODE{l}=NP;
end
% Boundary flags are added saved to NODE data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary flags for nodes are added to CELL data structure
for r=1:M
    P=CELL{r};
    for g=1:3
        ND=NODE{P{6+g}};
        P{9+g}=ND{2};
    end
    CELL{r}=P;
end
%% Boundary flags for nodes are added to CELL data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary flags for faces are added to CELL data structure
for r=1:M
    P=CELL{r};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 1
    FC=FACE{P{16}}; % Changing
    if FC{2}~=0 % The current face has no neighbor, so it is on boundary
        if P{10}==P{11} % The current face is in the middle of one type of boundary condition, changing
            P{19}=P{10}; % changing
        else % The current face is the joint of two type of boundary condition
            if (abs(P{10})/70)>=1 && (abs(P{11})/70)<1 % Changing
                P{19}=P{11}; % Changing
            elseif (abs(P{10})/70)<1 && (abs(P{11})/70)>=1 % Changing
                P{19}=P{10}; % Changing
            else
                error('The nodal boundary flags in the vacinity of corner nodes are incorrect!');
            end
        end
    else % The current face has neighbor, so it is interior
        P{19}=0; % Changing
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 2
    FC=FACE{P{17}}; % Changing
    if FC{2}~=0 % The current face has no neighbor, so it is on boundary
        if P{11}==P{12} % The current face is in the middle of one type of boundary condition, changing
            P{20}=P{11}; % changing
        else % The current face is the joint of two type of boundary condition
            if (abs(P{11})/70)>=1 && (abs(P{12})/70)<1 % Changing
                P{20}=P{12}; % Changing
            elseif (abs(P{11})/70)<1 && (abs(P{12})/70)>=1 % Changing
                P{20}=P{11}; % Changing
            else
                error('The nodal boundary flags in the vacinity of corner nodes are incorrect!');
            end
        end
    else % The current face has neighbor, so it is interior
        P{20}=0; % Changing
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 3
    FC=FACE{P{18}}; % Changing
    if FC{2}~=0 % The current face has no neighbor, so it is on boundary
        if P{12}==P{10} % The current face is in the middle of one type of boundary condition, changing
            P{21}=P{12}; % changing
        else % The current face is the joint of two type of boundary condition
            if (abs(P{12})/70)>=1 && (abs(P{10})/70)<1 % Changing
                P{21}=P{10}; % Changing
            elseif (abs(P{12})/70)<1 && (abs(P{10})/70)>=1 % Changing
                P{21}=P{12}; % Changing
            else
                error('The nodal boundary flags in the vacinity of corner nodes are incorrect!');
            end
        end
    else % The current face has neighbor, so it is interior
        P{21}=0; % Changing
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% face 3
    CELL{r}=P;
end

%% Boundary flags for faces are added to FACE data structure
for o=1:O
    FC=FACE{o};
    if FC{2}~=0
        neigh_cell_1=FC{12};
        neigh_cell_2=FC{13};
        if neigh_cell_1(1,1)==0 && neigh_cell_2(1,1)~=0
            neigh_cell=neigh_cell_2;
        elseif neigh_cell_1(1,1)~=0 && neigh_cell_2(1,1)==0
            neigh_cell=neigh_cell_1;
        else
            error('The neoghbor cell info for faces on the boundary contain false data!');
        end
        if neigh_cell(1,2)~=1 && neigh_cell(1,2)~=2 && neigh_cell(1,2)~=3
            error('The number could only be 1, 2 or 3!');
        end
        P=CELL{neigh_cell(1,1)};
        FC{2}=P{neigh_cell(1,2)+18};
    end
    FACE{o}=FC;
end

%% Find the normal stencil of corner nodes with in order to apply 
% periodic boundary conditions and extrapolation
% 1. One side is periodic, the other side is not
% ND{16}=[Exterior periodic node #;On-boundary periodic node #;Interior neighbor node #; Interior interior neighbor node #;..]
% 2. Both sides are periodic
% ND{16}=[Exterior periodic node #;On-boundary periodic node #;Interior neighbor node #; Interior interior neighbor node #;..;,
%        Exterior periodic node #;On-boundary periodic node #;Interior neighbor node #; Interior interior neighbor node #;..;]
% 3. Both sides are not periodic or One side is periodic, the other side is Fully Developed
% ND{16}=Empty
% Add more element in the column vector for larger stencil
for l=1:N
    NP=NODE{l};
    if (abs(NP{2})/70)>=1 % All corner nodes, including immersed corner nodes
        if abs(NP{2})==70 % Periodic + Periodic
            if NP{2}>0
                NS=NP{18};
                NU=NODE{NS(1)};
                ND=NODE{NS(2)};
                if NU{2}==1 && ND{2}==1
                    ;
                else
                    error('Check the boundary numerical flag on two sides of current corner node!');
                end
                for i=1:2
                    if i==1
                        Ncu=NU{3};
                        Ncd=NP{3};
                    else
                        Ncu=NP{3};
                        Ncd=ND{3};
                    end
                    L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
                    n_x=-(Ncd(2)-Ncu(2))/L;
                    n_y=(Ncd(1)-Ncu(1))/L;
                    n=[-n_x;-n_y];
                    if FM==0
                        NC=NP{3}+h*n;
                    elseif FM==1
                        NC=NP{3}+[dx*n(1,1);dy*n(2,1)];
                    else
                        error('Mesh flag is invalid!');
                    end
                    
                    XC=NC(1);
                    YC=NC(2);
                    if (XC>X2 || XC<X1) || (YC>Y2 || YC<Y1)
                        error('Reverse the normal vector!');
                    end
                    if abs(dot([1;0],n))>0 && single(e+abs(dot([0;1],n)))==single(e) % n is horizontal
                        N_temp=N_L;
                        if FM==0
                            spc=h;
                        elseif FM==1
                            spc=dx;
                        else
                            error('Mesh flag is invalid!');
                        end
                    elseif abs(dot([0;1],n))>0 && single(e+abs(dot([1;0],n)))==single(e) % n is vertical
                        N_temp=N_H;
                        if FM==0
                            spc=h;
                        elseif FM==1
                            spc=dy;
                        else
                            error('Mesh flag is invalid!');
                        end
                    else
                        error('The vector is not vertical or horizontal!');
                    end
                    %%%% Exterior neighbor
                    NC=NP{3}+(N_temp-2)*spc*n;
                    for k1=N_I:N
                        NT=NODE{k1};
                        NL=NC-NT{3};
                        XC=NL(1);
                        YC=NL(2);
                        if single(e+XC)==single(e) && single(e+YC)==single(e);
                            break;
                        end
                    end
                    if k1==N
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                    %%%% On-boundary neighbor
                    NC=NP{3}+(N_temp-1)*spc*n;
                    for k2=N_I:N
                        NT=NODE{k2};
                        NL=NC-NT{3};
                        XC=NL(1);
                        YC=NL(2);
                        if single(e+XC)==single(e) && single(e+YC)==single(e);
                            break;
                        end
                    end
                    if k2==N
                        NC=NP{3};
                        if (single(e+NC(1))==single(e+X2) && single(e+NC(2))==single(e+Y2)) || ((single(e+NC(1))==single(e+X1) && single(e+NC(2))==single(e+Y1)))
                            ;
                        else
                            error('On-boubdary neighbor for node on right wall not found!');
                        end
                    end
                    %%%% Interior neighbor
                    NC=NP{3}+spc*n;
                    for k3=1:N
                        NT=NODE{k3};
                        NL=NC-NT{3};
                        XC=NL(1);
                        YC=NL(2);
                        if single(e+XC)==single(e) && single(e+YC)==single(e);
                            break;
                        end
                    end
                    if k3==N
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                    %%%% Interior interior neighbor
                    NC=NP{3}+2*spc*n;
                    for k4=N_I:N
                        NT=NODE{k4};
                        NL=NC-NT{3};
                        XC=NL(1);
                        YC=NL(2);
                        if single(e+XC)==single(e) && single(e+YC)==single(e);
                            break;
                        end
                    end
                    if k4==N
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                    if i==1
                        NST(:,1)=[k1;k2;k3;k4];
                    else
                        NST(:,2)=[k1;k2;k3;k4];
                    end
                    NP{16}=NST;
                end
            else
                error('Periodic boundary is not available for immersed boundaries!');
            end
        elseif abs(NP{2})==71 || abs(NP{2})==72 % Periodic + Velocity or Velocity + Periodic || Periodic + Density or Density + Periodic 
            if NP{2}>0
                NS=NP{18};
                NU=NODE{NS(1)};
                ND=NODE{NS(2)};
                if NU{2}==1 && ND{2}~=1
                    Ncu=NU{3};
                    Ncd=NP{3};
                elseif NU{2}~=1 && ND{2}==1
                    Ncu=NP{3};
                    Ncd=ND{3};
                else
                    error('Check the boundary numerical flag on two sides of current corner node!');
                end
                L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
                n_x=-(Ncd(2)-Ncu(2))/L;
                n_y=(Ncd(1)-Ncu(1))/L;
                n=[-n_x;-n_y];
                if FM==0
                    NC=NP{3}+h*n;
                elseif FM==1
                    NC=NP{3}+[dx*n(1,1);dy*n(2,1)];
                else
                    error('Mesh flag is invalid!');
                end
                XC=NC(1);
                YC=NC(2);
                if (XC>X2 || XC<X1) || (YC>Y2 || YC<Y1)
                    error('Reverse the normal vector!');
                end
                if abs(dot([1;0],n))>0 && single(e+abs(dot([0;1],n)))==single(e) % n is horizontal
                    N_temp=N_L;
                    if FM==0
                        spc=h;
                    elseif FM==1
                        spc=dx;
                    else
                        error('Mesh flag is invalid!');
                    end
                elseif abs(dot([0;1],n))>0 && single(e+abs(dot([1;0],n)))==single(e) % n is vertical
                    N_temp=N_H;
                    if FM==0
                        spc=h;
                    elseif FM==1
                        spc=dy;
                    else
                        error('Mesh flag is invalid!');
                    end
                else
                    error('The vector is not vertical or horizontal!');
                end
                %%%% Exterior neighbor
                NC=NP{3}+(N_temp-2)*spc*n;
                for k1=N_I:N
                    NT=NODE{k1};
                    NL=NC-NT{3};
                    XC=NL(1);
                    YC=NL(2);
                    if single(e+XC)==single(e) && single(e+YC)==single(e);
                        break;
                    end
                end
                if k1==N
                    error('On-boubdary neighbor for node on right wall not found!');
                end
                %%%% On-boundary neighbor
                NC=NP{3}+(N_temp-1)*spc*n;
                for k2=N_I:N
                    NT=NODE{k2};
                    NL=NC-NT{3};
                    XC=NL(1);
                    YC=NL(2);
                    if single(e+XC)==single(e) && single(e+YC)==single(e);
                        break;
                    end
                end
                if k2==N
                    NC=NP{3};
                    if (single(e+NC(1))==single(e+X2) && single(e+NC(2))==single(e+Y2)) || (single(e+NC(1))==single(e+X1) && single(e+NC(2))==single(e+Y1))
                        ;
                    else
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                end
                %%%% Interior neighbor
                NC=NP{3}+spc*n;
                for k3=N_I:N
                    NT=NODE{k3};
                    NL=NC-NT{3};
                    XC=NL(1);
                    YC=NL(2);
                    if single(e+XC)==single(e) && single(e+YC)==single(e);
                        break;
                    end
                end
                if k3==N
                    error('On-boubdary neighbor for node on right wall not found!');
                end
                %%%% Interior interior neighbor
                NC=NP{3}+2*spc*n;
                for k4=N_I:N
                    NT=NODE{k4};
                    NL=NC-NT{3};
                    XC=NL(1);
                    YC=NL(2);
                    if single(e+XC)==single(e) && single(e+YC)==single(e);
                        break;
                    end
                end
                if k4==N
                    error('On-boubdary neighbor for node on right wall not found!');
                end
                %%%%%%% Stencil vector
                NP{16}=[k1;k2;k3;k4];
            else
                error('Periodic boundary is not available for immersed boundaries!');
            end
        elseif abs(NP{2})==73 % Velocity + Velocity
            ;
        elseif abs(NP{2})==74 % Density + Density
            ;
        elseif abs(NP{2})==75 % Velocity + Density or Density + Velocity  
            ;
        elseif abs(NP{2})==76 % Fully Developed + others or others
            ;
        elseif abs(NP{2})==77 % Fully Developed + Fully Developed
            ;
        else
            error('The numerical flag for current corner node is invalid or unavailable!');
        end
    else
       ;
    end
    NODE{l}=NP;
end
%% Find the periodic star structure of the nodes on boudaries in order to apply 
% % periodic boundary conditions
% % For any mesh as long as boundary nodes are lined up
% % ND{8} --- The total number of triangles connected to current periodic node
% % ND{9} --- The vector of order number of each connected triangle for
% % periodic boundary node (Empty for interior nodes)
% % ND{10} --- The vector of distance from the centroid of each connected
% % triangle to periodic boundary node (Empty for interior nodes)
% % ND{11} --- The sum of reverse distance from each star triangle to
% % periodic boundary node (Empty for interior nodes)
for l=1:N
    NP=NODE{l};
    if (abs(NP{2})/70)>=1 % All corner nodes, including inmmersed corner nodes
        if abs(NP{2})==70 % Periodic + Periodic
            if NP{2}>0
                NS=NP{16};
                NS1=NS(:,1);
                NS2=NS(:,2);
                for i=1:2
                    NC=0;
                    NDC=0;
                    if i==1
                        ND=NODE{NS1(2)};
                    else
                        ND=NODE{NS2(2)};
                    end
                    %%%% Total number of connected cells
                    if i==1
                        NP{8}=NP{4}+ND{4};
                    else
                        NP{8}=NP{8}+ND{4};
                    end
                    %%%% Vector of order # of connected cells
                    if i==1
                        NC=NP{5};
                    else
                        NC=NP{9};
                    end
                    L=length(ND{5});
                    NC(end+1:end+L)=ND{5};
                    NP{9}=NC;
                    if NP{8}~=length(NP{9})
                        error('The length of vector of order # of connected cell is wrong!');
                    end
                    %%%% Vector of distance
                    if i==1
                        NDC=NP{6};
                    else
                        NDC=NP{10};
                    end
                    L=length(ND{6});
                    NDC(end+1:end+L)=ND{6};
                    NP{10}=NDC;
                    if NP{8}~=length(NP{10})
                        error('The length of vector of distance is wrong!');
                    end
                end
                ND1=NODE{NS1(2)};
                ND2=NODE{NS2(2)};
                NSS1=ND1{16};
                NSS2=ND2{16};
                NSS11=NSS1(:,1);
                NSS22=NSS2(:,2);
                if NSS11(2)==NSS22(2) && NSS11(2)~=NP{1}
                    ;
                else
                    error('Did not find the diagonally facing corner node!');
                end
                ND=NODE{NSS11(2)};
                NC=0;
                NDC=0;
                %%%% Total number of connected cells
                NP{8}=NP{8}+ND{4};
                %%%% Vector of order # of connected cells
                NC=NP{9};
                L=length(ND{5});
                NC(end+1:end+L)=ND{5};
                NP{9}=NC;
                if NP{8}~=length(NP{9})
                    error('The length of vector of order # of connected cell is wrong!');
                end
                %%%% Vector of distance
                NDC=NP{10};
                L=length(ND{6});
                NDC(end+1:end+L)=ND{6};
                NP{10}=NDC;
                if NP{8}~=length(NP{10})
                    error('The length of vector of distance is wrong!');
                end
                %%%% Sum of reverse distance
                D_sum=0;
                for g=1:NP{8};
                    D_sum=D_sum+1/NDC(g);
                end;
                NP{11}=D_sum;
            else
                error('Periodic boundary is not available for immersed boundaries!');
            end
        elseif abs(NP{2})==71 || abs(NP{2})==72 % Periodic + Velocity or Velocity + Periodic || Periodic + Density or Density + Periodic 
            if NP{2}>0
                NC=0;
                NDC=0;
                NS=NP{16};
                ND=NODE{NS(2)};
                %%%% Total number of connected cells
                NP{8}=NP{4}+ND{4};
                %%%% Vector of order # of connected cells
                NC=NP{5};
                L=length(ND{5});
                NC(end+1:end+L)=ND{5};
                NP{9}=NC;
                if NP{8}~=length(NP{9})
                    error('The length of vector of order # of connected cell is wrong!');
                end
                %%%% Vector of distance
                NDC=NP{6};
                L=length(ND{6});
                NDC(end+1:end+L)=ND{6};
                NP{10}=NDC;
                if NP{8}~=length(NP{10})
                    error('The length of vector of distance is wrong!');
                end
                %%%% Sum of reverse distance
                D_sum=0;
                for g=1:NP{8};
                    D_sum=D_sum+1/NDC(g);
                end;
                NP{11}=D_sum;
            else
                error('Periodic boundary is not available for immersed boundaries!');
            end
        elseif abs(NP{2})==73 % Velocity + Velocity
            ;
        elseif abs(NP{2})==74 % Density + Density
            ;
        elseif abs(NP{2})==75 % Velocity + Density or Density + Velocity
            ;
        elseif abs(NP{2})==76 % Fully Developed + others or others
            ;
        elseif abs(NP{2})==77 % Fully Developed + Fully Developed
            ;
        else
            error('The numerical flag for current corner node is invalid or unavailable!');
        end
    else
       ;
    end
    NODE{l}=NP;
end
% %%%% Find the vector of pdf order # of upwind and bounceback of current
% %%%% node on boundary
% for l=1:N
%     ND=NODE{l};
%     if ND{2}==0 % Interior nodes
%         ;
%     elseif (abs(ND{2})/70)>=1 % All corner nodes
%         if ND{2}==70 % Periodic + Periodic
%             ;
%         elseif ND{2}==71 || ND{2}==72 % Periodic + Velocity or Velocity + Periodic || Periodic + Density or Density + Periodic 
%             Nng=ND{18};
%             NDU=NODE{Nng(1)};
%             NDD=NODE{Nng(2)};
%             if NDU{2}==1 && NDD{2}~=1
%                 Ncu=ND{3};
%                 Ncd=NDD{3};
%             elseif NDU{2}~=1 && NDD{2}==1
%                 Ncu=NDU{3};
%                 Ncd=ND{3};
%             else
%                 error('The current corner node should have been the joint of periodic and non-periodic BC!');
%             end
%             L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
%             n_x=-(Ncd(2)-Ncu(2))/L;
%             n_y=(Ncd(1)-Ncu(1))/L;
%             
%             UD1=zeros(1,q1);
%             for k=1:q1
%                 if dot([n_x,n_y],V1(:,k))>=0
%                     UD1(k)=k;
%                 else
%                     UD1(k)=0;
%                 end
%             end
%             for k=1:q1
%                 if UD1(k)==0
%                     for q=1:q1
%                         if single(e+V1(1,k)+V1(1,q))==single(e) && single(e+V1(2,k)+V1(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD1(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{19}=UD1;
%             
%             UD2=zeros(1,q2);
%             for k=1:q2
%                 if dot([n_x,n_y],V2(:,k))>=0
%                     UD2(k)=k;
%                 else
%                     UD2(k)=0;
%                 end
%             end
%             for k=1:q2
%                 if UD2(k)==0
%                     for q=1:q2
%                         if single(e+V2(1,k)+V2(1,q))==single(e) && single(e+V2(2,k)+V2(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD2(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{20}=UD2;
%         elseif abs(ND{2})==73 || abs(ND{2})==74 || abs(ND{2})==75 % Velocity + Velocity, Density + Density, Velocity + Density or Density + Velocity
%             Nng=ND{18};
%             NDU=NODE{Nng(1)};
%             NDD=NODE{Nng(2)};
%             Ncu=NDU{3};
%             Ncd=NDD{3};
%             L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
%             n_x=-(Ncd(2)-Ncu(2))/L;
%             n_y=(Ncd(1)-Ncu(1))/L;
%             
%             UD1=zeros(1,q1);
%             for k=1:q1
%                 if dot([n_x,n_y],V1(:,k))>=0
%                     UD1(k)=k;
%                 else
%                     UD1(k)=0;
%                 end
%             end
%             for k=1:q1
%                 if UD1(k)==0
%                     for q=1:q1
%                         if single(e+V1(1,k)+V1(1,q))==single(e) && single(e+V1(2,k)+V1(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD1(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{19}=UD1;
%             
%             UD2=zeros(1,q2);
%             for k=1:q2
%                 if dot([n_x,n_y],V2(:,k))>=0
%                     UD2(k)=k;
%                 else
%                     UD2(k)=0;
%                 end
%             end
%             for k=1:q2
%                 if UD2(k)==0
%                     for q=1:q2
%                         if single(e+V2(1,k)+V2(1,q))==single(e) && single(e+V2(2,k)+V2(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD2(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{20}=UD2;
%         elseif abs(ND{2})==76 % Fully Developed + others or others
%             ;
%         elseif abs(ND{2})==77 % Fully Developed + Fully Developed
%             ;
%         else
%             error('The numerical flag for current corner node is invalid or unavailable!');
%         end
%     else % Other boundary nodes
%         if (abs(ND{2})/70)==1
%             if ND{2}>0
%                 ;
%             else
%                 error('Immersed boundary cannot be periodic!');
%             end
%         else
%             Nng=ND{18};
%             NDU=NODE{Nng(1)};
%             NDD=NODE{Nng(2)};
%             Ncu=NDU{3};
%             Ncd=NDD{3};
%             L=sqrt((Ncu(1)-Ncd(1))^2+(Ncu(2)-Ncd(2))^2);
%             n_x=-(Ncd(2)-Ncu(2))/L;
%             n_y=(Ncd(1)-Ncu(1))/L;
%             
%             UD1=zeros(1,q1);
%             for k=1:q1
%                 if dot([n_x,n_y],V1(:,k))>=0
%                     UD1(k)=k;
%                 else
%                     UD1(k)=0;
%                 end
%             end
%             for k=1:q1
%                 if UD1(k)==0
%                     for q=1:q1
%                         if single(e+V1(1,k)+V1(1,q))==single(e) && single(e+V1(2,k)+V1(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD1(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{19}=UD1;
%             
%             UD2=zeros(1,q2);
%             for k=1:q2
%                 if dot([n_x,n_y],V2(:,k))>=0
%                     UD2(k)=k;
%                 else
%                     UD2(k)=0;
%                 end
%             end
%             for k=1:q2
%                 if UD2(k)==0
%                     for q=1:q2
%                         if single(e+V2(1,k)+V2(1,q))==single(e) && single(e+V2(2,k)+V2(2,q))==single(e)  % Two vectors are opposite
%                             break;
%                         end
%                     end
%                     UD2(k)=q;
%                 else
%                     ;
%                 end
%             end
%             ND{20}=UD2;
%         end
%     end
%     NODE{l}=ND;
% end