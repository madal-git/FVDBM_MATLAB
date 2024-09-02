function [CELL,NODE,FACE]=stencil(CELL,NODE,FACE,V1,V2,N_I,N_I_N,N_L,N_H,h,dx,dy,X,Y,X1,X2,Y1,Y2,FM)
% function [NODE,FACE]=stencil(CELL,NODE,FACE,V1,V2,N_I,N_I_N,N_L,N_H,h,dx,dy,X,Y,X1,X2,Y1,Y2,FM)
% fills all information that is stencil related into NODE and FACE.
% CELL is the cell data structure.
% NODE is the node data structure.
% FACE is the face data structure.
% V1 is the first lattice structure
% V2 is the second lattice structure
% X1, X2, Y1 and Y2 are bounds of mesh domain
% FM is flag for mesh type. FM=0----IRT mesh; FM=1----Random mesh

%% Data to be filled in NODE

% ND=NODE{r}

% ND{16}   The column vector of order # of the node stencil for outer plat boundaries(Empty for interior nodes)
%          , except for the corner nodes, which will be fill in bc.m. ********Empty for Triangle mesh**********

% ND{17}   The column vector of order # of the stencil for corner nodes along the diagonal direction (Empty for non-corner nodes).[near node #; further node #;..]********Empty for Triangle mesh**********
% ND{18}   The clockwise upstream and downstream boundary neighbor node of current boundary node(Empty for interior nodes). [upstream node #; downstream node #;..]
%% Data to be filled in FACE

% FC=FACE{r}

% *****************************************************************Before FACE v.11.24.2015*****************************************************************
% FC{16}   The assorted info of stencil points of current face based on V1
% FC{17}   The assorted info of stencil points of current face based on V2
% FC{18}   The coordinates of stencil points of current face based on the face normal
%          a 2_by_3 matrix = [coordinates of SP in downwind cell,coordinates of SP in upwind cell,coordinates of SP in further upwind cell]
% FC{19}   The Zone ID of stencil points of current face based on the face normal
%          a row vector = [Zone ID of SP in downwind cell, Zone ID of SP in upwind cell, Zone ID of SP in further upwind cell]

% Each item in FC{16} and FC{17} is a S=cell(25,1) and each element in S is a row vector listed as follows:
% {S{1}  downwind cell # at each direction of V1 or V2;
%  S{2}  upwind cell # at each direction of V1 or V2;
%  S{3}  further upwind cell # at each direction of V1 or V2;
%  S{4}  The x coordinates of downwind cell at each direction of V1 or V2;
%  S{5}  The y coordinates of downwind cell at each direction of V1 or V2;
%  S{6}  The x coordinates of upwind cell at each direction of V1 or V2;
%  S{7}  The y coordinates of upwind cell at each direction of V1 or V2;
%  S{8}  The x coordinates of further upwind cell at each direction of V1 or V2;
%  S{9}  The y coordinates of further upwind cell at each direction of V1 or V2;

%  S{10}  distance from downwind cell to upwind cell at each direction of V1 or V2;
%  S{11}  distance from upwind cell to further upwind cell at each direction of V1 or V2;
%  S{12}  One over the distance from downwind cell to upwind cell at each direction of V1 or V2;
%  S{13}  One over the distance from upwind cell to further upwind cell at each direction of V1 or V2;

%  S{14}  spatial ratio for current face,R1=(downwind cell size + upwind cell size)/upwind cell size,      at each direction of V1 or V2;
%  S{15}  spatial ratio for extrapolation, R2=distance from downwind cell to upwind cell/ distance from  upwind cell to further upwind cell, at each direction of V1 or V2;
%  S{16}  One over the spatial ratio for current face,R1=(downwind cell size + upwind cell size)/upwind cell size, at each direction of V1 or V2;
%  S{17}  One over the spatial ratio for extrapolation, R2=distance from downwind cell to upwind cell/ distance from  upwind cell to further upwind cell, at each direction of V1 or V2;

%  S{18} the left boundary node # on edge for extrapolation, at each direction of V1 or V2;
%  S{19} the right boundary node # on edge for extrapolation, at each direction of V1 or V2;
%  S{20} spatial ratio for on edge,R3=distance from left node to intercept node/distance from left node to right node], at each direction of V1 or V2}

% FC{20} The three points that circle the stencil point a 4_by_3 matrix = [three points for SP in downwind cell, three points for SP in upwind cell, three points for SP in further upwind cell]
% For Each column vector,
% [ S1=The total number of nodes in the three points:0-No node;1-one node;2-two nodes
%   S2=The order number of first point. If S1=0-order number of centroid; Else-order number of node
%   S3=The order number of second point. If S1<=1-order number of centroid; Else-order number of node
%   S4=The order number of third point. Always order of centroid, since S1<=2
% ]
% If S1=0, and S2=S3=S4, means the stencil is located at a centroid


%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{21}  The assorted info of 1st upwind of current face on boundary for V1
% FC{22}  The assorted info of 1st upwind of current face on boundary for V2

% FC{23}  The assorted info of 2nd upwind of current face on boundary for V1
% FC{24}  The assorted info of 2nd upwind of current face on boundary for V2
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% *****************************************************************Before FACE v.11.24.2015*****************************************************************

% *****************************************************************After FACE v.11.24.2015*****************************************************************
% FC{16} The coordinates of stencil points of current face based on the face normal a 2_by_4 matrix = [coordinates of SP in further downwind cell, coordinates of SP in downwind cell, coordinates of SP in upwind cell, coordinates of SP in further upwind cell]
% FC{17} The Zone ID of stencil points of current face based on the face normal
% a row vector = [Zone ID of SP in further downwind cell, Zone ID of SP in downwind cell, Zone ID of SP in upwind cell, Zone ID of SP in further upwind cell]
% FC{18} The three points that circle the stencil point a 4_by_4 matrix = [three points for SP in further downwind cell, three points for SP in downwind cell, three points for SP in upwind cell, three points for SP in further upwind cell]
% For Each column vector,
% [ S1=The total number of boundary nodes in the three points:0-No node;1-one node;2-two nodes
%  S2=The order number of first point. If S1=0-order number of centroid; Else-order number of first boundary node
%  S3=The order number of second point. If S1<=1-order number of centroid; Else-order number of second boundary node
%  S4=The order number of third point. Always order of centroid, since S1<=2
% ]
% If S1=0, and S2=S3=S4, means the stencil is located at a centroid

%%%%%%%%%%%%%%%%%%%%%%%%%%% Lattice-dependent Info%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{19}   The assorted info of stencil points of current face based on V1
% FC{20}   The assorted info of stencil points of current face based on V2
%% Each of them is a S=cell(30,1) and each element in S is a row vectors follows:

% {S{1}  downwind cell # at each direction of V1 or V2;
%  S{2}  upwind cell # at each direction of V1 or V2;
%  S{3}  further upwind cell # at each direction of V1 or V2;

%  S{4}  The x coordinates of stencil point in downwind cell at each direction of V1 or V2;
%  S{5}  The y coordinates of stencil point in downwind cell at each direction of V1 or V2;
%  S{6}  The x coordinates of stencil point in upwind cell at each direction of V1 or V2;
%  S{7}  The y coordinates of stencil point in upwind cell at each direction of V1 or V2;
%  S{8}  The x coordinates of stencil point in further upwind cell at each direction of V1 or V2;
%  S{9}  The y coordinates of stencil point in further upwind cell at each direction of V1 or V2;

%  S{10} The zone ID of stencil point in downwind cell at each direction of V1 or V2;
%  S{11} The zone ID of stencil point in upwind cell at each direction of V1 or V2;
%  S{12} The zone ID of stencil point in further upwind cell at each direction of V1 or V2;

%  S{13}  distance from the SP in downwind cell to the SP in upwind cell at each direction of V1 or V2;
%  S{14}  distance from the SP in upwind cell to the SP in further upwind cell at each direction of V1 or V2;
%  S{15}  One over S{13};
%  S{16}  One over S{14};

%  S{17}  spatial ratio for current face,R1=(downwind cell size + upwind cell size)/upwind cell size,      at each direction of V1 or V2;
%  S{18}  spatial ratio for extrapolation, R2=distance from the SP in downwind cell to the SP in upwind cell/ distance from the SP in upwind cell to the SP in further upwind cell, at each direction of V1 or V2;
%  S{19}  One over S{17};
%  S{20}  One over S{18};

%  S{21} the left boundary node # on edge for extrapolation, at each direction of V1 or V2;
%  S{22} the right boundary node # on edge for extrapolation, at each direction of V1 or V2;
%  S{23} spatial ratio for on edge,R3=distance from left node to intercept node/distance from left node to right node], at each direction of V1 or V2}
% FC{21}  The upwind identifier for 2nd-order mapping at the stencil points for V1
% a q1_by_4 matrix (q1 the number of velocities in V1) = [upwind identifier for SP in further downwind cell, upwind identifier for SP in downwind cell, upwind identifier for SP in upwind cell, upwind identifier for SP in further upwind cell]
% For Each column vector, each entry is either 0 or 1
% 0---The SP is located upwind of the centroid, use 2nd-order mapping for SP
% 1---The SP is located downwind of the centroid, use the value at centroid for SP
% FC{22}  The upwind identifier for 2nd-order mapping at the stencil points for V2
% a q2_by_4 matrix (q2 the number of velocities in V2) = [upwind identifier for SP in further downwind cell, upwind identifier for SP in downwind cell, upwind identifier for SP in upwind cell, upwind identifier for SP in further upwind cell]
%%%%%%%%%%%%%%%%%%%%%%%%%%% Lattice-dependent Info%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC{23}  The coordinates of stencil points of current face based on the face normal if one or many of the stencil points is out of boundaries
% FC{24}  The Zone ID of stencil points of current face based on the face normal if one or many of the stencil points is out of boundaries
% FC{25}  The three points that circle the stencil point if one or many of the stencil points is out of boundaries
% FC{26}  The assorted info of stencil points of current face based on V1 if one or many of the stencil points is out of boundaries
% FC{27}  The assorted info of stencil points of current face based on V2 if one or many of the stencil points is out of boundaries
% FC{28}  The upwind identifier for 2nd-order mapping at the stencil points for V1 if one or many of the stencil points is out of boundaries
% FC{29}  The upwind identifier for 2nd-order mapping at the stencil points for V2 if one or many of the stencil points is out of boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%% For Periodic boundaries ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%

% *****************************************************************After FACE v.11.24.2015*****************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=length(CELL);
N=length(NODE);
O=length(FACE);
q1=length(V1(1,:));
q2=length(V2(1,:));
K=30;
P=CELL{1};
e=ceil(10*sqrt(M*P{6})); %%%% Calculate the reference size
e1=e*100;
Tol=1e-6; %%%% Tolerence for finding the intercept point on the edge for 2nd-order upwind scheme
% The last cell in CELL
E=CELL{M};
Eb=0;
for s=1:3
    if E{19+s-1}~=0
        Eb=Eb+1;
    end
end
if Eb~=0
    msg=['The last cell has ', num2str(Eb), ' faces on the boundary!'];
    disp(msg);
end
G_counter=0;
BL=[X1;Y1];
BR=[X2;Y1];
TL=[X1;Y2];
TR=[X2;Y2];
%% Fill in the data in NODE
%% ND{16} filling
ND=NODE{N-1};
if isempty(ND{16})
    for l=1:N
        if l<=N_I
            ;
        elseif l==N_I+N_L-1 || l==N_I+N_L-1+N_H-1 || l==N_I+N_L-1+N_H-2+N_L || l==N % To be done in bc.m
            ;
        else
            ND=NODE{l};
            if l<=N_I+N_L-2 %%%% Top
                if FM==0 %%%% IRT mesh
                    %%%% Exterior neighbor
                    for k1=1:N
                        if single(e+X(l))==single(e+X(k1)) && single(e+Y(l))==single(e+Y(k1)+(N_H-2)*h);
                            break;
                        end
                    end
                    if k1==N
                        error('Exterior neighbor for node on top wall not found!');
                    end
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+X(l))==single(e+X(k2)) && single(e+Y(l))==single(e+Y(k2)+(N_H-1)*h);
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on top wall not found!');
                    end
                    %%%% Interior neighbor
                    for k3=1:N
                        if single(e+X(l))==single(e+X(k3)) && single(e+Y(l))==single(e+Y(k3)+h);
                            break;
                        end
                    end
                    if k3==N
                        error('Interior neighbor for node on top wall not found!');
                    end
                    %%%% Interior interior neighbor
                    for k4=1:N
                        if single(e+X(l))==single(e+X(k4)) && single(e+Y(l))==single(e+Y(k4)+2*h);
                            break;
                        end
                    end
                    if k4==N
                        error('Interior neighbor for node on top wall not found!');
                    end
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                elseif FM==1 % Random mesh
                    %%%% Exterior neighbor
                    k1=0;
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+X(l))==single(e+X(k2)) && single(e+Y(l))==single(e+Y(k2)+(N_H-1)*dy);
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on top wall not found!');
                    end
                    %%%% Interior neighbor
                    k3=0;
                    %%%% Interior interior neighbor
                    k4=0;
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                else
                    error('Mesh marker FM has to be 0 or 1!');
                end
            elseif l<=N_I+N_L-1+N_H-2 %%%% Right
                if FM==0 %%%% IRT mesh
                    %%%% Exterior neighbor
                    for k1=1:N
                        if single(e+Y(l))==single(e+Y(k1)) && single(e+X(l))==single(e+X(k1)+(N_L-2)*h);
                            break;
                        end
                    end
                    if k1==N
                        error('Exterior neighbor for node on right wall not found!');
                    end
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+Y(l))==single(e+Y(k2)) && single(e+X(l))==single(e+X(k2)+(N_L-1)*h);
                            break;
                        end
                    end
                    if k2==N
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                    %%%% Interior neighbor
                    for k3=1:N
                        if single(e+Y(l))==single(e+Y(k3)) && single(e+X(l))==single(e+X(k3)+h);
                            break;
                        end
                    end
                    if k3==N
                        error('Interior neighbor for node on right wall not found!');
                    end
                    %%%% Interior interior neighbor
                    for k4=1:N
                        if single(e+Y(l))==single(e+Y(k4)) && single(e+X(l))==single(e+X(k4)+2*h);
                            break;
                        end
                    end
                    if k4==N
                        error('Interior neighbor for node on right wall not found!');
                    end
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                elseif FM==1 % Random mesh
                    %%%% Exterior neighbor
                    k1=0;
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+Y(l))==single(e+Y(k2)) && single(e+X(l))==single(e+X(k2)+(N_L-1)*dx);
                            break;
                        end
                    end
                    if k2==N
                        error('On-boubdary neighbor for node on right wall not found!');
                    end
                    %%%% Interior neighbor
                    k3=0;
                    %%%% Interior interior neighbor
                    k4=0;
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                else
                    error('Mesh marker FM has to be 0 or 1!');
                end
            elseif l<=N_I+N_L-1+N_H-2+N_L-1 %%%% Bottom
                if FM==0 %%%% IRT mesh
                    %%%% Exterior neighbor
                    for k1=1:N
                        if single(e+X(l))==single(e+X(k1)) && single(e+Y(l)+(N_H-2)*h)==single(e+Y(k1));
                            break;
                        end
                    end
                    if k1==N
                        error('Exterior neighbor for node on bottom wall not found!');
                    end
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+X(l))==single(e+X(k2)) && single(e+Y(l)+(N_H-1)*h)==single(e+Y(k2));
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on bottom wall not found!');
                    end
                    %%%% Interior neighbor
                    for k3=1:N
                        if single(e+X(l))==single(e+X(k3)) && single(e+Y(l)+h)==single(e+Y(k3));
                            break;
                        end
                    end
                    if k3==N
                        error('Interior neighbor for node on bottom wall not found!');
                    end
                    %%%% Interior interior neighbor
                    for k4=1:N
                        if single(e+X(l))==single(e+X(k4)) && single(e+Y(l)+2*h)==single(e+Y(k4));
                            break;
                        end
                    end
                    if k4==N
                        error('Interior neighbor for node on bottom wall not found!');
                    end
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                elseif FM==1 % Random mesh
                    %%%% Exterior neighbor
                    k1=0;
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+X(l))==single(e+X(k2)) && single(e+Y(l)+(N_H-1)*dy)==single(e+Y(k2));
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on bottom wall not found!');
                    end
                    %%%% Interior neighbor
                    k3=0;
                    %%%% Interior interior neighbor
                    k4=0;
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                else
                    error('Mesh marker FM has to be 0 or 1!');
                end
            else %%%% Left
                if FM==0 %%%% IRT mesh
                    %%%% Exterior neighbor
                    for k1=1:N
                        if single(e+Y(l))==single(e+Y(k1)) && single(e+X(l)+(N_L-2)*h)==single(e+X(k1));
                            break;
                        end
                    end
                    if k1==N
                        error('Exterior neighbor for node on left wall not found!');
                    end
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+Y(l))==single(e+Y(k2)) && single(e+X(l)+(N_L-1)*h)==single(e+X(k2));
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on left wall not found!');
                    end
                    %%%% Interior neighbor
                    for k3=1:N
                        if single(e+Y(l))==single(e+Y(k3)) && single(e+X(l)+h)==single(e+X(k3));
                            break;
                        end
                    end
                    if k3==N
                        error('Exterior neighbor for node on left wall not found!');
                    end
                    %%%% Interior interior neighbor
                    for k4=1:N
                        if single(e+Y(l))==single(e+Y(k4)) && single(e+X(l)+2*h)==single(e+X(k4));
                            break;
                        end
                    end
                    if k4==N
                        error('Exterior neighbor for node on left wall not found!');
                    end
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                elseif FM==1 % Random mesh
                    %%%% Exterior neighbor
                    k1=0;
                    %%%% On-boundary neighbor
                    for k2=N_I:N
                        if single(e+Y(l))==single(e+Y(k2)) && single(e+X(l)+(N_L-1)*dx)==single(e+X(k2));
                            break;
                        end
                    end
                    if k2==N
                        error('On-boundary neighbor for node on left wall not found!');
                    end
                    %%%% Interior neighbor
                    k3=0;
                    %%%% Interior interior neighbor
                    k4=0;
                    %%%%%%% Stencil vector
                    ND{16,1}=[k1;k2;k3;k4];
                else
                    error('Mesh marker FM has to be 0 or 1!');
                end
            end
            NODE{l}=ND;
        end
    end
end

%% ND{17} Filling, only for IRT mesh
if FM==0
    ND=NODE{N};
    if isempty(ND{17})
        for l=1:N
            ND=NODE{l};
            if l==N_I+N_L-1 %%%% Top right corner
                %%% First node(nearer)
                for k=1:N
                    if single(e+X(l))==single(e+X(k)+h/2) && single(e+Y(l))==single(e+Y(k)+h/2);
                        break;
                    end
                end
                if k==N
                    error('First diagonal node for top-right corner node is not found!');
                end
                %%% Second node(further)
                for r=1:N
                    if single(e+X(l))==single(e+X(r)+h) && single(e+Y(l))==single(e+Y(r)+h);
                        break;
                    end
                end
                if r==N
                    error('Second diagonal node for top-right corner node is not found!');
                end
                ND{17,1}=[k;r];
                NODE{l}=ND;
            elseif l==N_I+N_L-1+N_H-1 %%%% Bottom right  corner
                %%% First node(nearer)
                for k=1:N
                    if single(e+X(l))==single(e+X(k)+h/2) && single(e+Y(l))==single(e+Y(k)-h/2);
                        break;
                    end
                end
                if k==N
                    error('First diagonal node for bottom-right corner node is not found!');
                end
                %%% Second node(further)
                for r=1:N
                    if single(e+X(l))==single(e+X(r)+h) && single(e+Y(l))==single(e+Y(r)-h);
                        break;
                    end
                end
                if r==N
                    error('Second diagonal node for bottom-right corner node is not found!');
                end
                ND{17,1}=[k;r];
                NODE{l}=ND;
            elseif l==N_I+N_L-1+N_H-2+N_L %%%% Bottom left  corner
                %%% First node(nearer)
                for k=1:N
                    if single(e+X(l))==single(e+X(k)-h/2) && single(e+Y(l))==single(e+Y(k)-h/2);
                        break;
                    end
                end
                if k==N
                    error('First diagonal node for bottom-left corner node is not found!');
                end
                %%% Second node(further)
                for r=1:N
                    if single(e+X(l))==single(e+X(r)-h) && single(e+Y(l))==single(e+Y(r)-h);
                        break;
                    end
                end
                if r==N
                    error('Second diagonal node for bottom-left corner node is not found!');
                end
                ND{17,1}=[k;r];
                NODE{l}=ND;
            elseif l==N %%%% Top left corner
                %%% First node(nearer)
                for k=1:N
                    if single(e+X(l))==single(e+X(k)-h/2) && single(e+Y(l))==single(e+Y(k)+h/2);
                        break;
                    end
                end
                if k==N
                    error('First diagonal node for top-left corner node is not found!');
                end
                %%% Second node(further)
                for r=1:N
                    if single(e+X(l))==single(e+X(r)-h) && single(e+Y(l))==single(e+Y(r)+h);
                        break;
                    end
                end
                if r==N
                    error('Second diagonal node for top-left corner node is not found!');
                end
                ND{17,1}=[k;r];
                NODE{l}=ND;
            else
                ;
            end
            NODE{l}=ND;
        end
    end
end
%% ND{18} filling
% %%%%%  Clockwiselly the upstream and downstream boundary neighbor node for current outer boundary node
% %%%%%  Or, counter-clockwiselly the upstream and downstream boundary neighbor node for current inner boundary node
% %%%%%  A column vector with first row for upstream node and second row the downstream node
% %%%%%  (Empty for interior nodes)
ND=NODE{N};
if isempty(ND{18})
    for l=1:N
        if l<=N_I
            if l<=N_I_N  % Immersed boundary nodes, current version is only for one hole
                ND=NODE{l};
                if l==N_I_N %%%% The last immersed boundary node
                    Um=l-1;
                    Dm=1;
                    ND{18,1}=[Um;Dm];
                elseif l==1 %%%% The first immersed boundary node
                    Um=N_I_N;
                    Dm=l+1;
                    ND{18,1}=[Um;Dm];
                else
                    Um=l-1;
                    Dm=l+1;
                    ND{18,1}=[Um;Dm];
                end
                NODE{l}=ND;
            else
                ;
            end
        else  % Outer boundary nodes
            ND=NODE{l};
            if l==N_I+N_L-1 %%%% Top right corner
                Um=l-1;
                Dm=l+1;
                ND{18,1}=[Um;Dm];
            elseif l==N_I+N_L-1+N_H-1 %%%% Right bottom corner
                Um=l-1;
                Dm=l+1;
                ND{18,1}=[Um;Dm];
            elseif l==N_I+N_L-1+N_H-2+N_L %%%% Left bottom corner
                Um=l-1;
                Dm=l+1;
                ND{18,1}=[Um;Dm];
            elseif l==N %%%% Top left corner
                Um=l-1;
                Dm=N_I+1;
                ND{18,1}=[Um;Dm];
            elseif l==N_I+1 %%%% The first boundary node
                Um=N;
                Dm=l+1;
                ND{18,1}=[Um;Dm];
            else
                Um=l-1;
                Dm=l+1;
                ND{18,1}=[Um;Dm];
            end
            NODE{l}=ND;
        end
    end
end


%% FACE info
FC=FACE{1};
if isempty(FC{16}) || isempty(FC{17})
    % Fill in FC{16}, {17}, {18}
    for l=1:O
        FC=FACE{l};
        neigh_up=FC{12};
        neigh_down=FC{13};
        C_b=FC{7}; % The base point coordinates of the stencil
        
        %% FC{16}, {17}, {18}
        fc16=[0,neigh_down(1,1),neigh_up(1,1),0];
        fc17=zeros(2,4);
        fc18=zeros(3,4);
        if fc16(1,2)==0 && fc16(1,3)~=0 % The face is on boundary, the downwind cell is out of boundary
            fc16(1,1)=0; % If the downwind cell is out of boundary, the further downwind cell is also out of boundary
            cell_up=CELL{fc16(1,3)};
            Nd1=NODE{cell_up{7}};
            Nd2=NODE{cell_up{8}};
            Nd3=NODE{cell_up{9}};
            %% Find the coordinate of upwind stencil point, since the
            % further downwind and downwind stencil points are defaulted to
            % be zero
            C_j_u=norm_joint(C_b,FC{4}',cell_up{5});
            if in_triangle(C_j_u,Nd1{3},Nd2{3},Nd3{3})
                fc17(:,3)=C_j_u;
            else
                error('The generated point is not within the current cell!');
            end
            %% Find the coordinate of further upwind stencil point
            nd_group_uu=[cell_up{7},cell_up{8},cell_up{9}];
            nd_opp_face_uu=setxor([FC{8},FC{9}],nd_group_uu); %find the node opposing the face
            if length(nd_opp_face_uu)~=1
                error('The node opposing the face is not found!');
            end
            Nd_opp_face_uu=NODE{nd_opp_face_uu};
            cell_group_uu=setxor(cell_up{1},Nd_opp_face_uu{5});
            if (Nd_opp_face_uu{4}-length(cell_group_uu))~=1
                error('The upwind cell should be excluded!');
            end
            c_j_uu=zeros(2,Nd_opp_face_uu{4}-1);
            uu_found_counter=0;
            uu_found_cell=0;
            c_j_uu_found=zeros(2,1);
            for i=1:Nd_opp_face_uu{4}-1
                P=CELL{cell_group_uu(i)};
                Nd1=NODE{P{7}};
                Nd2=NODE{P{8}};
                Nd3=NODE{P{9}};
                c_j_uu(:,i)=norm_joint(C_b,FC{4}',P{5});
                if in_triangle(c_j_uu(:,i),Nd1{3},Nd2{3},Nd3{3})
                    uu_found_counter=uu_found_counter+1;
                    uu_found_cell(uu_found_counter)=cell_group_uu(i);
                    c_j_uu_found(:,uu_found_counter)=c_j_uu(:,i);
                end
            end
            if uu_found_counter==0
                error('The further upwind cell is not found!');
            elseif uu_found_counter==1
                if uu_found_cell==0
                    error('Logic error!');
                end
                fc16(1,4)=uu_found_cell;
                fc17(:,4)=c_j_uu_found;
            else % More than one further upwind cell is found
                Dist=zeros(1,uu_found_counter);
                for i=1:uu_found_counter
                    Dist(i)=dis(C_b,c_j_uu_found(:,i));
                end
                Dist_min=min(Dist);
                for i=1:uu_found_counter
                    if single(e+Dist(i))==single(e+Dist_min)
                        break;
                    end
                end
                if i==uu_found_counter
                    if single(e+Dist(i))~=single(e+Dist_min)
                        error('Logic error!');
                    end
                end
                fc16(1,4)=uu_found_cell(i);
                fc17(:,4)=c_j_uu_found(:,i);
            end
            % Check
            cell_found=CELL{fc16(1,4)};
            Nd1=NODE{cell_found{7}};
            Nd2=NODE{cell_found{8}};
            Nd3=NODE{cell_found{9}};
            if ~in_triangle(fc17(:,4),Nd1{3},Nd2{3},Nd3{3})
                error('The further upwind stencil point is not within the found cell!');
            end
            if single(e+dis(fc17(:,4),fc17(:,3))+dis(fc17(:,3),C_b))~=single(e+dis(fc17(:,4),C_b))
                error('The stencil points are not aligned!');
            end
            %% fill in FC{17} and {18} for downwind and further downwind stencil points
            fc17(:,2)=C_b+FC{4}'*dis(C_b,fc17(:,3));
            fc18(:,2)=[FC{8};FC{9};0.5];
            fc17(:,1)=fc17(:,2)+FC{4}'*dis(C_b,fc17(:,2));
            fc18(:,1)=[FC{8};FC{9};0.5];
        elseif fc16(1,3)==0 && fc16(1,2)~=0 % The face is on boundary, the upwind cell is out of boundary
            fc16(1,4)=0; % If the upwind cell is out of boundary, the further upwind cell is also out of boundary
            cell_down=CELL{fc16(1,2)};
            Nd1=NODE{cell_down{7}};
            Nd2=NODE{cell_down{8}};
            Nd3=NODE{cell_down{9}};
            %% Find the coordinate of downwind stencil point, since the
            % further upwind and upwind stencil points are defaulted to
            % be zero
            C_j_d=norm_joint(C_b,FC{4}',cell_down{5});
            if in_triangle(C_j_d,Nd1{3},Nd2{3},Nd3{3})
                fc17(:,2)=C_j_d;
            else
                error('The generated point is not within the current cell!');
            end
            %% Find the coordinate of further downwind stencil point
            nd_group_dd=[cell_down{7},cell_down{8},cell_down{9}];
            nd_opp_face_dd=setxor([FC{8},FC{9}],nd_group_dd); %find the node opposing the face
            if length(nd_opp_face_dd)~=1
                error('The node opposing the face is not found!');
            end
            Nd_opp_face_dd=NODE{nd_opp_face_dd};
            cell_group_dd=setxor(cell_down{1},Nd_opp_face_dd{5});
            if (Nd_opp_face_dd{4}-length(cell_group_dd))~=1
                error('The downwind cell should be excluded!');
            end
            c_j_dd=zeros(2,Nd_opp_face_dd{4}-1);
            dd_found_counter=0;
            dd_found_cell=0;
            c_j_dd_found=zeros(2,1);
            for i=1:Nd_opp_face_dd{4}-1
                P=CELL{cell_group_dd(i)};
                Nd1=NODE{P{7}};
                Nd2=NODE{P{8}};
                Nd3=NODE{P{9}};
                c_j_dd(:,i)=norm_joint(C_b,FC{4}',P{5});
                if in_triangle(c_j_dd(:,i),Nd1{3},Nd2{3},Nd3{3})
                    dd_found_counter=dd_found_counter+1;
                    dd_found_cell(dd_found_counter)=cell_group_dd(i);
                    c_j_dd_found(:,dd_found_counter)=c_j_dd(:,i);
                end
            end
            if dd_found_counter==0
                error('The further upwind cell is not found!');
            elseif dd_found_counter==1
                if dd_found_cell==0
                    error('Logic error!');
                end
                fc16(1,1)=dd_found_cell;
                fc17(:,1)=c_j_dd_found;
            else % More than one further upwind cell is found
                Dist=zeros(1,dd_found_counter);
                for i=1:dd_found_counter
                    Dist(i)=dis(C_b,c_j_dd_found(:,i));
                end
                Dist_min=min(Dist);
                for i=1:dd_found_counter
                    if single(e+Dist(i))==single(e+Dist_min)
                        break;
                    end
                end
                if i==dd_found_counter
                    if single(e+Dist(i))~=single(e+Dist_min)
                        error('Logic error!');
                    end
                end
                fc16(1,1)=dd_found_cell(i);
                fc17(:,1)=c_j_dd_found(:,i);
            end
            % Check
            cell_found=CELL{fc16(1,1)};
            Nd1=NODE{cell_found{7}};
            Nd2=NODE{cell_found{8}};
            Nd3=NODE{cell_found{9}};
            if ~in_triangle(fc17(:,1),Nd1{3},Nd2{3},Nd3{3})
                error('The further downwind stencil point is not within the found cell!');
            end
            if single(e+dis(fc17(:,1),fc17(:,2))+dis(fc17(:,2),C_b))~=single(e+dis(fc17(:,1),C_b))
                error('The stencil points are not aligned!');
            end
            %% fill in FC{17} and {18} for upwind and further upwind stencil points
            fc17(:,3)=C_b-FC{4}'*dis(C_b,fc17(:,2));
            fc18(:,3)=[FC{8};FC{9};0.5];
            fc17(:,4)=fc17(:,3)-FC{4}'*dis(C_b,fc17(:,3));
            fc18(:,4)=[FC{8};FC{9};0.5];
        elseif fc16(1,2)~=0 && fc16(1,3)~=0
            cell_up=CELL{fc16(1,3)};
            Nd1_up=NODE{cell_up{7}};
            Nd2_up=NODE{cell_up{8}};
            Nd3_up=NODE{cell_up{9}};
            cell_down=CELL{fc16(1,2)};
            Nd1_down=NODE{cell_down{7}};
            Nd2_down=NODE{cell_down{8}};
            Nd3_down=NODE{cell_down{9}};
            %% Find the coordinate of upwind stencil point
            C_j_u=norm_joint(C_b,FC{4}',cell_up{5});
            if in_triangle(C_j_u,Nd1_up{3},Nd2_up{3},Nd3_up{3})
                fc17(:,3)=C_j_u;
            else
                error('The generated point is not within the current cell!');
            end
            %% Find the coordinate of downwind stencil point
            C_j_d=norm_joint(C_b,FC{4}',cell_down{5});
            if in_triangle(C_j_d,Nd1_down{3},Nd2_down{3},Nd3_down{3})
                fc17(:,2)=C_j_d;
            else
                error('The generated point is not within the current cell!');
            end
            %% Check
            if single(e+dis(fc17(:,2),C_b)+dis(fc17(:,3),C_b))~=single(e+dis(fc17(:,2),fc17(:,3)))
                error('The stencil points are not aligned!');
            end
            %% Find the coordinate of further upwind stencil point
            nd_group_uu=[cell_up{7},cell_up{8},cell_up{9}];
            nd_opp_face_uu=setxor([FC{8},FC{9}],nd_group_uu); %find the node opposing the face
            if length(nd_opp_face_uu)~=1
                error('The node opposing the face is not found!');
            end
            Nd_opp_face_uu=NODE{nd_opp_face_uu};
            cell_group_uu=setxor(cell_up{1},Nd_opp_face_uu{5});
            if (Nd_opp_face_uu{4}-length(cell_group_uu))~=1
                error('The upwind cell should be excluded!');
            end
            c_j_uu=zeros(2,Nd_opp_face_uu{4}-1);
            uu_found_counter=0;
            uu_found_cell=0;
            c_j_uu_found=zeros(2,1);
            for i=1:Nd_opp_face_uu{4}-1
                P=CELL{cell_group_uu(i)};
                Nd1=NODE{P{7}};
                Nd2=NODE{P{8}};
                Nd3=NODE{P{9}};
                c_j_uu(:,i)=norm_joint(C_b,FC{4}',P{5});
                if in_triangle(c_j_uu(:,i),Nd1{3},Nd2{3},Nd3{3})
                    uu_found_counter=uu_found_counter+1;
                    uu_found_cell(uu_found_counter)=cell_group_uu(i);
                    c_j_uu_found(:,uu_found_counter)=c_j_uu(:,i);
                end
            end
            if uu_found_counter==0
                the_other_two_face=setxor(FC{1},[cell_up{16},cell_up{17},cell_up{18}]);
                if length(the_other_two_face)~=2
                    error('Logic Error!');
                end
                FC1=FACE{the_other_two_face(1)};
                FC2=FACE{the_other_two_face(2)};
                if intersect([FC1{8},FC1{9}],[FC2{8},FC2{9}])~=nd_opp_face_uu
                    error('Logic Error!');
                end
                [Bool_found_1, C_i_1] = intercept(C_b,FC{4}',FC1{5},FC1{6});
                [Bool_found_2, C_i_2] = intercept(C_b,FC{4}',FC2{5},FC2{6});
                if Nd_opp_face_uu{2}==0
                    error('The opposing node must be on boundary!');
                else
                    if Bool_found_1==1 && Bool_found_2==0
                        if FC1{2}==0 % The intercept point must be extended to another boundary face that shares the opposing node
                            nd_neigh=Nd_opp_face_uu{18};
                            Nd_up=NODE{nd_neigh(1)};
                            Nd_down=NODE{nd_neigh(2)};
                            [Bool_found_backup_1, C_i_backup_1]=intercept(C_b,FC{4}',Nd_up{3},Nd_opp_face_uu{3});
                            [Bool_found_backup_2, C_i_backup_2]=intercept(C_b,FC{4}',Nd_down{3},Nd_opp_face_uu{3});
                            if Bool_found_backup_1==1 && Bool_found_backup_2==0
                                fc17(:,4)=C_i_backup_1-FC{4}'*dis(C_i_backup_1,fc17(:,3));
                                fc18(:,4)=[nd_opp_face_uu;nd_neigh(1);dis(C_i_backup_1,Nd_opp_face_uu{3})/dis(Nd_up{3},Nd_opp_face_uu{3})];
                            elseif Bool_found_backup_1==0 && Bool_found_backup_2==1
                                fc17(:,4)=C_i_backup_2-FC{4}'*dis(C_i_backup_2,fc17(:,3));
                                fc18(:,4)=[nd_opp_face_uu;nd_neigh(2);dis(C_i_backup_2,Nd_opp_face_uu{3})/dis(Nd_down{3},Nd_opp_face_uu{3})];
                            else
                                if Bool_found_backup_1==0 && Bool_found_backup_2==0
                                    error('Check mesh!');
                                else
                                    error('Logic Error!');
                                end
                            end
                        else % the intercept point in on FC1
                            fc17(:,4)=C_i_1-FC{4}'*dis(C_i_1,fc17(:,3));
                            fc18(:,4)=[FC1{8};FC1{9};dis(C_i_1,FC1{5})/FC1{3}];
                        end
                    elseif Bool_found_1==0 && Bool_found_2==1
                        if FC2{2}==0 % The intercept point must be extended to another boundary face that shares the opposing node
                            nd_neigh=Nd_opp_face_uu{18};
                            Nd_up=NODE{nd_neigh(1)};
                            Nd_down=NODE{nd_neigh(2)};
                            [Bool_found_backup_1, C_i_backup_1]=intercept(C_b,FC{4}',Nd_up{3},Nd_opp_face_uu{3});
                            [Bool_found_backup_2, C_i_backup_2]=intercept(C_b,FC{4}',Nd_down{3},Nd_opp_face_uu{3});
                            if Bool_found_backup_1==1 && Bool_found_backup_2==0
                                fc17(:,4)=C_i_backup_1-FC{4}'*dis(C_i_backup_1,fc17(:,3));
                                fc18(:,4)=[nd_opp_face_uu;nd_neigh(1);dis(C_i_backup_1,Nd_opp_face_uu{3})/dis(Nd_up{3},Nd_opp_face_uu{3})];
                            elseif Bool_found_backup_1==0 && Bool_found_backup_2==1
                                fc17(:,4)=C_i_backup_2-FC{4}'*dis(C_i_backup_2,fc17(:,3));
                                fc18(:,4)=[nd_opp_face_uu;nd_neigh(2);dis(C_i_backup_2,Nd_opp_face_uu{3})/dis(Nd_down{3},Nd_opp_face_uu{3})];
                            else
                                if Bool_found_backup_1==0 && Bool_found_backup_2==0
                                    error('Check mesh!');
                                else
                                    error('Logic Error!');
                                end
                            end
                        else % the intercept point is on FC2
                            fc17(:,4)=C_i_2-FC{4}'*dis(C_i_2,fc17(:,3));
                            fc18(:,4)=[FC2{8};FC2{9};dis(C_i_2,FC2{5})/FC2{3}];
                        end
                    else
                        if Bool_found_1==1 && Bool_found_2==1
                            if single(e+sum(C_i_1-C_i_2))==single(e) % The stencil is right through the opposing node
                                fc17(:,4)=Nd_opp_face_uu{3}-FC{4}'*dis(Nd_opp_face_uu{3},fc17(:,3));
                                fc18(:,4)=[nd_opp_face_uu;nd_opp_face_uu;0.5];
                            else
                                error('Logic error!');
                            end
                        else
                            error('Logic error!');
                        end
                    end
                end
                fc16(1,4)=0;
            elseif uu_found_counter==1
                if uu_found_cell==0
                    error('Logic error!');
                end
                fc16(1,4)=uu_found_cell;
                fc17(:,4)=c_j_uu_found;
            else % More than one further upwind cell is found
                Dist=zeros(1,uu_found_counter);
                for i=1:uu_found_counter
                    Dist(i)=dis(C_b,c_j_uu_found(:,i));
                end
                Dist_min=min(Dist);
                for i=1:uu_found_counter
                    if single(e+Dist(i))==single(e+Dist_min)
                        break;
                    end
                end
                if i==uu_found_counter
                    if single(e+Dist(i))~=single(e+Dist_min)
                        error('Logic error!');
                    end
                end
                fc16(1,4)=uu_found_cell(i);
                fc17(:,4)=c_j_uu_found(:,i);
            end
            % Check
            if fc16(1,4)~=0
                cell_found=CELL{fc16(1,4)};
                Nd1=NODE{cell_found{7}};
                Nd2=NODE{cell_found{8}};
                Nd3=NODE{cell_found{9}};
                if ~in_triangle(fc17(:,4),Nd1{3},Nd2{3},Nd3{3})
                    error('The further upwind stencil point is not within the found cell!');
                end
                if single(e+dis(fc17(:,4),fc17(:,3))+dis(fc17(:,3),C_b))~=single(e+dis(fc17(:,4),C_b))
                    error('The stencil points are not aligned!');
                end
            end
            %% Find the coordinate of further downwind stencil point
            nd_group_dd=[cell_down{7},cell_down{8},cell_down{9}];
            nd_opp_face_dd=setxor([FC{8},FC{9}],nd_group_dd); %find the node opposing the face
            if length(nd_opp_face_dd)~=1
                error('The node opposing the face is not found!');
            end
            Nd_opp_face_dd=NODE{nd_opp_face_dd};
            cell_group_dd=setxor(cell_down{1},Nd_opp_face_dd{5});
            if (Nd_opp_face_dd{4}-length(cell_group_dd))~=1
                error('The downwind cell should be excluded!');
            end
            c_j_dd=zeros(2,Nd_opp_face_dd{4}-1);
            dd_found_counter=0;
            dd_found_cell=0;
            c_j_dd_found=zeros(2,1);
            for i=1:Nd_opp_face_dd{4}-1
                P=CELL{cell_group_dd(i)};
                Nd1=NODE{P{7}};
                Nd2=NODE{P{8}};
                Nd3=NODE{P{9}};
                c_j_dd(:,i)=norm_joint(C_b,FC{4}',P{5});
                if in_triangle(c_j_dd(:,i),Nd1{3},Nd2{3},Nd3{3})
                    dd_found_counter=dd_found_counter+1;
                    dd_found_cell(dd_found_counter)=cell_group_dd(i);
                    c_j_dd_found(:,dd_found_counter)=c_j_dd(:,i);
                end
            end
            if dd_found_counter==0
                the_other_two_face=setxor(FC{1},[cell_down{16},cell_down{17},cell_down{18}]);
                if length(the_other_two_face)~=2
                    error('Logic Error!');
                end
                FC1=FACE{the_other_two_face(1)};
                FC2=FACE{the_other_two_face(2)};
                if intersect([FC1{8},FC1{9}],[FC2{8},FC2{9}])~=nd_opp_face_dd
                    error('Logic Error!');
                end
                
                [Bool_found_1, C_i_1] = intercept(C_b,FC{4}',FC1{5},FC1{6});
                [Bool_found_2, C_i_2] = intercept(C_b,FC{4}',FC2{5},FC2{6});
                
                if Nd_opp_face_dd{2}==0
                    FC{1}
                    Nd_opp_face_dd{1}
                    cell_group_dd
                    fc_1=FC{5};
                    fc_2=FC{6};
                        plot([fc_1(1),fc_2(1)],[fc_1(2),fc_2(2)],'red', 'linewidth',1);

                        hold on;
    pause(3)
                    
for r=1:Nd_opp_face_dd{4}-1;
    P=CELL{cell_group_dd(r)};
    plot(P{22},P{23},'black', 'linewidth',1);
    hold on;
    pause(0.2)
    plot(P{24},P{25},'black', 'linewidth',1);
    hold on;
    pause(0.2)
    plot(P{26},P{27},'black', 'linewidth',1);
    hold on;
    pause(0.2)
    CenD=P{5};
    
    plot(CenD(1),CenD(2),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'blue');
            pause(2)
        hold on
        c_c=norm_joint(C_b,FC{4}',P{5});
        plot(c_c(1),c_c(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'm');
                pause(2)
        hold on
end;
    axis equal tight;

    
                    error('The opposing node must be on boundary!');
                else
                    if Bool_found_1==1 && Bool_found_2==0
                        if FC1{2}==0
                            nd_neigh=Nd_opp_face_dd{18};
                            Nd_up=NODE{nd_neigh(1)};
                            Nd_down=NODE{nd_neigh(2)};
                            [Bool_found_backup_1, C_i_backup_1]=intercept(C_b,FC{4}',Nd_up{3},Nd_opp_face_dd{3});
                            [Bool_found_backup_2, C_i_backup_2]=intercept(C_b,FC{4}',Nd_down{3},Nd_opp_face_dd{3});
                            if Bool_found_backup_1==1 && Bool_found_backup_2==0
                                fc17(:,1)=C_i_backup_1+FC{4}'*dis(C_i_backup_1,fc17(:,2));
                                fc18(:,1)=[nd_opp_face_dd;nd_neigh(1);dis(C_i_backup_1,Nd_opp_face_dd{3})/dis(Nd_up{3},Nd_opp_face_dd{3})];
                            elseif Bool_found_backup_1==0 && Bool_found_backup_2==1
                                fc17(:,1)=C_i_backup_2+FC{4}'*dis(C_i_backup_2,fc17(:,2));
                                fc18(:,1)=[nd_opp_face_dd;nd_neigh(2);dis(C_i_backup_2,Nd_opp_face_dd{3})/dis(Nd_down{3},Nd_opp_face_dd{3})];
                            else
                                if Bool_found_backup_1==0 && Bool_found_backup_2==0
                                    error('Check mesh!');
                                else
                                    error('Logic Error!');
                                end
                            end
                        else
                            fc17(:,1)=C_i_1+FC{4}'*dis(C_i_1,fc17(:,2));
                            fc18(:,1)=[FC1{8};FC1{9};dis(C_i_1,FC1{5})/FC1{3}];
                        end
                    elseif Bool_found_1==0 && Bool_found_2==1
                        if FC2{2}==0
                            nd_neigh=Nd_opp_face_dd{18};
                            Nd_up=NODE{nd_neigh(1)};
                            Nd_down=NODE{nd_neigh(2)};
                            [Bool_found_backup_1, C_i_backup_1]=intercept(C_b,FC{4}',Nd_up{3},Nd_opp_face_dd{3});
                            [Bool_found_backup_2, C_i_backup_2]=intercept(C_b,FC{4}',Nd_down{3},Nd_opp_face_dd{3});
                            if Bool_found_backup_1==1 && Bool_found_backup_2==0
                                fc17(:,1)=C_i_backup_1+FC{4}'*dis(C_i_backup_1,fc17(:,2));
                                fc18(:,1)=[nd_opp_face_dd;nd_neigh(1);dis(C_i_backup_1,Nd_opp_face_dd{3})/dis(Nd_up{3},Nd_opp_face_dd{3})];
                            elseif Bool_found_backup_1==0 && Bool_found_backup_2==1
                                fc17(:,1)=C_i_backup_2+FC{4}'*dis(C_i_backup_2,fc17(:,2));
                                fc18(:,1)=[nd_opp_face_dd;nd_neigh(2);dis(C_i_backup_2,Nd_opp_face_dd{3})/dis(Nd_down{3},Nd_opp_face_dd{3})];
                            else
                                if Bool_found_backup_1==0 && Bool_found_backup_2==0
                                    error('Check mesh!');
                                else
                                    error('Logic Error!');
                                end
                            end
                        else
                            fc17(:,1)=C_i_2+FC{4}'*dis(C_i_2,fc17(:,2));
                            fc18(:,1)=[FC2{8};FC2{9};dis(C_i_2,FC2{5})/FC2{3}];
                        end
                    else
                        if Bool_found_1==1 && Bool_found_2==1
                            if single(e+sum(C_i_1-C_i_2))==single(e) % The stencil is right through the opposing node
                                fc17(:,1)=Nd_opp_face_dd{3}+FC{4}'*dis(Nd_opp_face_dd{3},fc17(:,2));
                                fc18(:,1)=[nd_opp_face_dd;nd_opp_face_dd;0.5];
                            else
                                error('Logic error!');
                            end
                        else
                            error('Logic error!');
                        end
                    end
                end
                fc16(1,1)=0;
            elseif dd_found_counter==1
                if dd_found_cell==0
                    error('Logic error!');
                end
                fc16(1,1)=dd_found_cell;
                fc17(:,1)=c_j_dd_found;
            else % More than one further upwind cell is found
                Dist=zeros(1,dd_found_counter);
                for i=1:dd_found_counter
                    Dist(i)=dis(C_b,c_j_dd_found(:,i));
                end
                Dist_min=min(Dist);
                for i=1:dd_found_counter
                    if single(e+Dist(i))==single(e+Dist_min)
                        break;
                    end
                end
                if i==dd_found_counter
                    if single(e+Dist(i))~=single(e+Dist_min)
                        error('Logic error!');
                    end
                end
                fc16(1,1)=dd_found_cell(i);
                fc17(:,1)=c_j_dd_found(:,i);
            end
            % Check
            if fc16(1,1)~=0
                cell_found=CELL{fc16(1,1)};
                Nd1=NODE{cell_found{7}};
                Nd2=NODE{cell_found{8}};
                Nd3=NODE{cell_found{9}};
                if ~in_triangle(fc17(:,1),Nd1{3},Nd2{3},Nd3{3})
                    error('The further upwind stencil point is not within the found cell!');
                end
                if single(e+dis(fc17(:,1),fc17(:,2))+dis(fc17(:,2),C_b))~=single(e+dis(fc17(:,1),C_b))
                    error('The stencil points are not aligned!');
                end
            end
        else
            error('The face must be either on boundary or interior!');
        end
        FC{16}=fc16;
        FC{17}=fc17;
        FC{18}=fc18;
        %% Final data filling
        FACE{l}=FC;
    end
    % Final check for FC{16}, {17} and {18}
    counter1=0;
    counter2=0;
    counter3=0;
    counter4=0;
    counter5=0;
    counter6=0;
    for l=1:O
        FC=FACE{l};
        C_b=FC{7};
        fc16=FC{16};
        fc17=FC{17};
        fc18=FC{18};
        if single(e+dis(fc17(:,1),fc17(:,4)))~=single(e+dis(fc17(:,1),fc17(:,2))+dis(fc17(:,2),C_b)+dis(C_b,fc17(:,3))+dis(fc17(:,3),fc17(:,4)))
            error('The stencil points are not aligned!');
        end
        
        if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose further downwind stencil point is out of boundary
            if FC{2}~=0
                error('This should be an interior face!');
            end
            counter1=counter1+1;
        end
        
        if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose further upwind stencil point is out of boundary
            if FC{2}~=0
                error('This should be an interior face!');
            end
            counter2=counter2+1;
        end
        
        if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Check the face whose both further downwind  and further upwind stencil points are out of boundary
            if FC{2}~=0
                error('This should be an interior face!');
            end
            counter3=counter3+1;
        end
        
        if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check the face whose all stencil points are within the boundary
            if FC{2}~=0
                error('This should be an interior face!');
            end
            counter4=counter4+1;
        end
        
        if (fc16(1,1)==0 && fc16(1,2)==0) && (fc16(1,3)~=0 && fc16(1,4)~=0) % Check boundary face whose all downwind side stencil points are out of boundary
            if FC{2}==0
                error('This should be an interior face!');
            end
            counter5=counter5+1;
        end
        
        if (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)==0 && fc16(1,4)==0) % Check boundary face whose all upwind side stencil points are out of boundary
            if FC{2}==0
                error('This should be an interior face!');
            end
            counter6=counter6+1;
        end
        for i=1:4
            if fc16(1,i)==0
                Nd1=NODE{fc18(1,i)};
                Nd2=NODE{fc18(2,i)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i=Nd1{3}+fc18(3,i)*L*v;
                if FC{2}==0
                    if Nd1{1}==Nd2{1}
                        if L~=0;
                            error('The intercept point for out-of-boundary stencil points are incorrect!');
                        end
                    else
                        n=C_i-C_b;
                        n=n/norm(n);
                        if single(e+abs(FC{4}*n))~=single(e+1)
                            error('The intercept point for out-of-boundary stencil points are incorrect!');
                        end
                    end
                else
                    if single(e+sum(C_i-C_b))~=single(e)
                        error('The intercept point for out-of-boundary stencil points are incorrect!');
                    end
                end
            end
        end
        
    end
    if (counter1+counter2+counter3+counter4+counter5+counter6)~=O
        error('Logic error!');
    end
    
    %% FC{19}
    e1=1e-7;
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc19=zeros(1,length(fc16));
        for i=1:length(fc16)
            if fc16(1,i)~=0
                cell_in=CELL{fc16(1,i)};
                Nd1=NODE{cell_in{7}};
                Nd2=NODE{cell_in{8}};
                Nd3=NODE{cell_in{9}};
                C_nd1=Nd1{3};
                C_nd2=Nd2{3};
                C_nd3=Nd3{3};
                C_c=cell_in{5};
                if dis(fc17(:,i),C_c)<e1
                    in_zone_four=1;
                else
                    in_zone_one=in_triangle(fc17(:,i),C_c,C_nd1,C_nd2);
                    in_zone_two=in_triangle(fc17(:,i),C_c,C_nd2,C_nd3);
                    in_zone_three=in_triangle(fc17(:,i),C_c,C_nd3,C_nd1);
                    in_zone_four=0;
                end
                % Check
                if (in_zone_one+in_zone_two+in_zone_three)~=1 && (in_zone_one+in_zone_two+in_zone_three)~=0
                    error('The stencil point could no be in multiple zones!');
                end
                % Further determination
                if in_zone_four
                    fc19(1,i)=4;
                else
                    if in_zone_one
                        fc19(1,i)=1; % Located in zone 1
                    elseif in_zone_two
                        fc19(1,i)=2; % Located in zone 2
                    elseif in_zone_three
                        fc19(1,i)=3; % Located in zone 3
                    else
                        on_edge_one=on_edge(fc17(:,i),C_c,C_nd1);
                        on_edge_two=on_edge(fc17(:,i),C_c,C_nd2);
                        on_edge_three=on_edge(fc17(:,i),C_c,C_nd3);
                        % Check
                        if (on_edge_one+on_edge_two+on_edge_three)~=1 && (on_edge_one+on_edge_two+on_edge_three)~=3
                            error('The stencil point could no be on multiple edges!');
                        end
                        if on_edge_one && (on_edge_two && on_edge_three)
                            % Check
                            if single(e+dis(C_c,fc17(:,i)))~=single(e)
                                error('Logic error!');
                            end
                            fc19(1,i)=4;
                        elseif on_edge_one
                            fc19(1,i)=1;
                        elseif on_edge_two
                            fc19(1,i)=2;
                        elseif on_edge_three
                            fc19(1,i)=3;
                        else
                            error('The zone for the stencil point is not found!');
                        end
                    end
                end
            end
        end
        FC{19}=fc19;
        %% Final data filling
        FACE{l}=FC;
    end
    % Check for FC{19}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc19=FC{19};
        for i=1:length(fc16)
            if fc16(1,i)~=0
                Cell_fix=CELL{fc16(1,i)};
                C_fix=Cell_fix{5};
                Nd1=NODE{Cell_fix{7}};
                Nd2=NODE{Cell_fix{8}};
                Nd3=NODE{Cell_fix{9}};
                if fc19(1,i)==1
                    C1=Nd1{3};
                    C2=Nd2{3};
                    if on_edge(fc17(:,i),C_fix,C1);
                        break;
                    end
                elseif fc19(1,i)==2
                    C1=Nd2{3};
                    C2=Nd3{3};
                    if on_edge(fc17(:,i),C_fix,C1);
                        break;
                    end
                elseif fc19(1,i)==3
                    C1=Nd3{3};
                    C2=Nd1{3};
                    if on_edge(fc17(:,i),C_fix,C1);
                        break;
                    end
                elseif fc19(1,i)==4
                    if dis(fc17(:,i),C_fix)>e1
                        error('The stencil point should be at the centroid!');
                    else
                        break;
                    end
                else
                    error('The points that enlose the stencil point cannot be all nodes!');
                end
                if ~in_triangle(fc17(:,i),C_fix,C1,C2)
                    error('Wrong zone ID!');
                end
            end
        end
    end
    
    
    %% FC{20}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc19=FC{19};
        fc20=zeros(4,length(fc16));
        for i=1:length(fc16)
            if fc16(1,i)~=0
                if fc19(1,i)==4
                    fc20(1,i)=0;
                    fc20(2,i)=fc16(1,i);
                    fc20(3,i)=fc16(1,i);
                    fc20(4,i)=fc16(1,i);
                else
                    Cell_fix=CELL{fc16(1,i)};
                    Face_target=FACE{Cell_fix{15+fc19(1,i)}};
                    if Face_target{10}~=0 && Face_target{11}~=0 % Both end nodes of target face are located on boundary
                        if Face_target{2}==0
                            error('The outer and inner boundaries are too close!');
                        end
                        fc20(1,i)=2;
                        fc20(2,i)=Face_target{8};
                        fc20(3,i)=Face_target{9};
                        fc20(4,i)=fc16(1,i);
                    elseif (Face_target{10}~=0 && Face_target{11}==0) || (Face_target{10}==0 && Face_target{11}~=0) % One of the two end nodes of target face is located on boundary
                        % Check
                        if Face_target{2}~=0
                            error('The current face should be interior!');
                        end
                        if Face_target{10}~=0 && Face_target{11}==0
                            ND=NODE{Face_target{8}};
                            C_nd=ND{3};
                            nd=Face_target{8};
                        elseif Face_target{10}==0 && Face_target{11}~=0
                            ND=NODE{Face_target{9}};
                            C_nd=ND{3};
                            nd=Face_target{9};
                        else
                            error('There must be one node on the boundary and the other is not!');
                        end
                        % Algorithm 1, always use the boundary node
                        %% Find all possible cells that could be used to form triangles
                        C_fixed=Cell_fix{5};
                        ND1=NODE{Face_target{8}};
                        ND2=NODE{Face_target{9}};
                        Cell_third_pool=union(ND1{5},ND2{5});
                        Cell_third_pool=setxor(fc16(1,i),Cell_third_pool);
                        L=length(Cell_third_pool);
                        a=0;
                        in_tri=0;
                        for k=1:L
                            Cell_third=CELL{Cell_third_pool(k)};
                            C_third_cell=Cell_third{5};
                            if in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
                                a=a+1;
                                in_tri(a)=k;
                            end
                        end
                        if a==0
                            error('Enclosing triangle is not found!');
                        end
                        %% Find the cell that has the shortest distance to the stencil point
                        Dis=zeros(1,a);
                        for k=1:a
                            Cell_third=CELL{Cell_third_pool(in_tri(k))};
                            C_third_cell=Cell_third{5};
                            Dis(1,k)=dis(fc17(:,i),C_third_cell);
                        end
                        Dis_min=min(Dis);
                        for k=1:a
                            if double(e+Dis_min)==double(e+Dis(1,k))
                                break;
                            end
                        end
                        if k==a
                            if Dis(1,k)~=Dis_min
                                error('Logic error!');
                            end
                        end
                        Third_cell_found=Cell_third_pool(in_tri(k));
                        % Check
                        Cell_third=CELL{Third_cell_found};
                        C_third_cell=Cell_third{5};
                        if ~in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
                            error('Enclosing triangle is not found!');
                        end
                        %%
                        fc20(1,i)=1;
                        fc20(2,i)=nd;
                        fc20(3,i)=fc16(1,i);
                        fc20(4,i)=Third_cell_found;
                        
                        
%                         %% Algorithm 2, try centroids, if no found, then use
%                         % the boundary node
%                         %% Find all possible paired centroids that could be used to form triangles
%                         C_fixed=Cell_fixed{5};
%                         ND1=NODE{Face_target{8}};
%                         ND2=NODE{Face_target{9}};
%                         Cell_union=union(ND1{5},ND2{5});
%                         Cell_union=setxor(fc16(1,i),Cell_union);
%                         L=length(Cell_union);
%                         Pair_cell_union=zeros(2,L*(L-1)/2);
%                         a=0;
%                         for k=1:length(Cell_union)
%                             if length(Cell_union)==2
%                                 a=a+1;
%                                 Pair_cell_union(:,a)=[Cell_union(1);Cell_union(2)];
%                                 break;
%                             else
%                                 for n=2:length(Cell_union)
%                                     a=a+1;
%                                     Pair_cell_union(:,a)=[Cell_union(1);Cell_union(n)];
%                                 end
%                                 Cell_union=setxor(Cell_union(1),Cell_union);
%                             end
%                         end
%                         % check
%                         if a~=L*(L-1)/2
%                             error('Logic error!');
%                         end
%                         %% Narrow down the previous possibility ruling out the triangles that don't circle the stencil point in the upstream cell
%                         b=0;
%                         in_tri_pair=0;
%                         for k=1:a
%                             Cell_pair_1=CELL{Pair_cell_union(1,k)};
%                             Cell_pair_2=CELL{Pair_cell_union(2,k)};
%                             C_1=Cell_pair_1{5};
%                             C_2=Cell_pair_2{5};
%                             if in_triangle(fc17(:,i),C_fixed,C_1,C_2)
%                                 b=b+1;
%                                 in_tri_pair(b)=k;
%                             end
%                         end
%                         if b==0
%                             %% Find all possible cells that could be used to form triangles
%                             ND1=NODE{Face_target{8}};
%                             ND2=NODE{Face_target{9}};
%                             Cell_third_pool=union(ND1{5},ND2{5});
%                             Cell_third_pool=setxor(fc16(1,i),Cell_third_pool);
%                             L=length(Cell_third_pool);
%                             a=0;
%                             in_tri=0;
%                             for k=1:L
%                                 Cell_third=CELL{Cell_third_pool(k)};
%                                 C_third_cell=Cell_third{5};
%                                 if in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
%                                     a=a+1;
%                                     in_tri(a)=k;
%                                 end
%                             end
%                             if a==0
%                                 error('Enclosing triangle is not found!');
%                             end
%                             %% Find the cell that has the shortest distance to the stencil point
%                             Dis=zeros(1,a);
%                             for k=1:a
%                                 Cell_third=CELL{Cell_third_pool(in_tri(k))};
%                                 C_third_cell=Cell_third{5};
%                                 Dis(1,k)=dis(fc17(:,i),C_third_cell);
%                             end
%                             Dis_min=min(Dis);
%                             for k=1:a
%                                 if single(e+Dis_min)==single(e+Dis(1,k))
%                                     break;
%                                 end
%                             end
%                             if k==a
%                                 if Dis(1,k)~=Dis_min
%                                     error('Logic error!');
%                                 end
%                             end
%                             Third_cell_found=Cell_third_pool(in_tri(k));
%                             % Check
%                             Cell_third=CELL{Third_cell_found};
%                             C_third_cell=Cell_third{5};
%                             if ~in_triangle(fc17(:,i),C_fixed,C_nd,C_third_cell)
%                                 error('Enclosing triangle is not found!');
%                             end
%                             %%
%                             fc20(1,i)=1;
%                             fc20(2,i)=nd;
%                             fc20(3,i)=fc16(1,i);
%                             fc20(4,i)=Third_cell_found;
%                         else
%                             %% Find the cell that has the shortest distance to the stencil point
%                             Dis=zeros(1,b);
%                             for k=1:b
%                                 Cell_pair_1=CELL{Pair_cell_union(1,in_tri_pair(k))};
%                                 Cell_pair_2=CELL{Pair_cell_union(2,in_tri_pair(k))};
%                                 Dis(1,k)=(dis(fc17(:,i),Cell_pair_1{5})+dis(fc17(:,i),Cell_pair_2{5}))/2;
%                             end
%                             Dis_min=min(Dis);
%                             for k=1:b
%                                 if single(e+Dis_min)==single(e+Dis(1,k))
%                                     break;
%                                 end
%                             end
%                             if k==b
%                                 if Dis(1,k)~=Dis_min
%                                     error('Logic error!');
%                                 end
%                             end
%                             %% Check
%                             Cell_pair_1=CELL{Pair_cell_union(1,in_tri_pair(k))};
%                             Cell_pair_2=CELL{Pair_cell_union(2,in_tri_pair(k))};
%                             C_1=Cell_pair_1{5};
%                             C_2=Cell_pair_2{5};
%                             if ~in_triangle(fc17(:,i),C_fixed,C_1,C_2)
%                                 error('Logic error!');
%                             end
%                             fc20(1,i)=0;
%                             fc20(2,i)=fc16(:,i);
%                             fc20(3,i)=Pair_cell_union(1,in_tri_pair(k));
%                             fc20(4,i)=Pair_cell_union(2,in_tri_pair(k));
%                         end
                    elseif Face_target{10}==0 && Face_target{11}==0 % Both end nodes of target face are interior nodes
                        C_fixed=Cell_fix{5};
                        %% Find all possible paired centroids that could be used to form triangles
                        ND1=NODE{Face_target{8}};
                        ND2=NODE{Face_target{9}};
                        Cell_union=union(ND1{5},ND2{5});
                        Cell_union=setxor(fc16(1,i),Cell_union);
                        L=length(Cell_union);
                        Pair_cell_union=zeros(2,L*(L-1)/2);
                        a=0;
                        for k=1:length(Cell_union)
                            if length(Cell_union)==2
                                a=a+1;
                                Pair_cell_union(:,a)=[Cell_union(1);Cell_union(2)];
                                break;
                            else
                                for n=2:length(Cell_union)
                                    a=a+1;
                                    Pair_cell_union(:,a)=[Cell_union(1);Cell_union(n)];
                                end
                                Cell_union=setxor(Cell_union(1),Cell_union);
                            end
                        end
                        % check
                        if a~=L*(L-1)/2
                            error('Logic error!');
                        end
                        %% Narrow down the previous possibility ruling out the triangles that don't circle the stencil point in the upstream cell
                        b=0;
                        in_tri_pair=0;
                        for k=1:a
                            Cell_pair_1=CELL{Pair_cell_union(1,k)};
                            Cell_pair_2=CELL{Pair_cell_union(2,k)};
                            C_1=Cell_pair_1{5};
                            C_2=Cell_pair_2{5};
                            if in_triangle(fc17(:,i),C_fixed,C_1,C_2)
                                b=b+1;
                                in_tri_pair(b)=k;
                            end
                        end
                        if b==0
                            error('No triangle is found that contains the stencil point!');
                        end
                        %% Find the cell that has the shortest distance to the stencil point
                        Dis=zeros(1,b);
                        for k=1:b
                            Cell_pair_1=CELL{Pair_cell_union(1,in_tri_pair(k))};
                            Cell_pair_2=CELL{Pair_cell_union(2,in_tri_pair(k))};
                            Dis(1,k)=(dis(fc17(:,i),Cell_pair_1{5})+dis(fc17(:,i),Cell_pair_2{5}))/2;
                        end
                        Dis_min=min(Dis);
                        for k=1:b
                            if double(e+Dis_min)==double(e+Dis(1,k))
                                break;
                            end
                        end
                        if k==b
                            if Dis(1,k)~=Dis_min
                                error('Logic error!');
                            end
                        end
                        %% Check
                        Cell_pair_1=CELL{Pair_cell_union(1,in_tri_pair(k))};
                        Cell_pair_2=CELL{Pair_cell_union(2,in_tri_pair(k))};
                        C_1=Cell_pair_1{5};
                        C_2=Cell_pair_2{5};
                        if ~in_triangle(fc17(:,i),C_fixed,C_1,C_2)
                            error('Logic error!');
                        end
                        fc20(1,i)=0;
                        fc20(2,i)=fc16(:,i);
                        fc20(3,i)=Pair_cell_union(1,in_tri_pair(k));
                        fc20(4,i)=Pair_cell_union(2,in_tri_pair(k));
                    else
                        error('Logic error!');
                    end
                end
            end
        end
        FC{20}=fc20;
        %% Final data filling
        FACE{l}=FC;
    end
    % Check for FC{20}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc20=FC{20};
        for i=1:length(fc16)
            if fc16(1,i)~=0
                if fc20(1,i)==0
                    if fc20(2,i)==fc20(3,i) && fc20(3,i)==fc20(4,i)
                        
                    else
                        Cell1=CELL{fc20(2,i)};
                        Cell2=CELL{fc20(3,i)};
                        Cell3=CELL{fc20(4,i)};
                        if ~in_triangle(fc17(:,i),Cell1{5},Cell2{5},Cell3{5})
                            error('Wrong enclosing points!');
                        end
                    end
                elseif fc20(1,i)==1
                    Node1=NODE{fc20(2,i)};
                    Cell2=CELL{fc20(3,i)};
                    Cell3=CELL{fc20(4,i)};
                    if ~in_triangle(fc17(:,i),Node1{3},Cell2{5},Cell3{5})
                        error('Wrong enclosing points!');
                    end
                elseif fc20(1,i)==2
                    Node1=NODE{fc20(2,i)};
                    Node2=NODE{fc20(3,i)};
                    Cell3=CELL{fc20(4,i)};
                    if ~in_triangle(fc17(:,i),Node1{3},Node2{3},Cell3{5})
                        error('Wrong enclosing points!');
                    end
                else
                    error('The points that enlose the stencil point cannot be all nodes!');
                end
            end
        end
    end
    
    %% FC{21}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc20=FC{20};
        fc21=zeros(6,length(fc16));
        for i=1:length(fc16)
            if fc16(1,i)~=0
                if fc20(1,i)==0
                    Cell1=CELL{fc20(2,i)};
                    Cell2=CELL{fc20(3,i)};
                    Cell3=CELL{fc20(4,i)};
                    fc21(1:2,i)=Cell1{5};
                    fc21(3:4,i)=Cell2{5};
                    fc21(5:6,i)=Cell3{5};
                elseif fc20(1,i)==1
                    Node1=NODE{fc20(2,i)};
                    Cell2=CELL{fc20(3,i)};
                    Cell3=CELL{fc20(4,i)};
                    fc21(1:2,i)=Node1{3};
                    fc21(3:4,i)=Cell2{5};
                    fc21(5:6,i)=Cell3{5};
                elseif fc20(1,i)==2
                    Node1=NODE{fc20(2,i)};
                    Node2=NODE{fc20(3,i)};
                    Cell3=CELL{fc20(4,i)};
                    fc21(1:2,i)=Node1{3};
                    fc21(3:4,i)=Node2{3};
                    fc21(5:6,i)=Cell3{5};
                else
                    error('Logic error');
                end
            end
        end
        FC{21}=fc21;
        FACE{l}=FC;
    end
    % Check FC{21}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc19=FC{19};
        fc21=FC{21};
        for i=1:length(fc16)
            if fc16(1,i)~=0
                C1=fc21(1:2,i);
                C2=fc21(3:4,i);
                C3=fc21(5:6,i);
                if fc19(1,i)==4
                    if single(e+dis(C1,fc17(:,i)))~=single(e) || (single(e+dis(C2,fc17(:,i)))~=single(e) || single(e+dis(C3,fc17(:,i)))~=single(e))
                        error('FC{21} contains false info!');
                    end
                else
                    if ~in_triangle(fc17(:,i),C1,C2,C3)
                        error('FC{21} contains false info!');
                    end
                end
            end
        end
    end
    
    %% FC{22}
    for l=1:O
        FC=FACE{l};
        C_b=FC{7};
        fc16=FC{16};
        fc17=FC{17};
        fc22=zeros(1,length(fc16));
        for i=1:length(fc16)
            fc22(1,i)=dis(C_b,fc17(:,i));
        end
        FC{22}=fc22;
        FACE{l}=FC;
    end
    % check FC{22}
    for l=1:O
        FC=FACE{l};
        fc17=FC{17};
        fc22=FC{22};
        if fc22(1,1)<fc22(1,2)
            error('The further downwind stencil point should be futher away from the face than the downwind stencil point!');
        end
        if fc22(1,4)<fc22(1,3)
            error('The further upwind stencil point should be futher away from the face than the upwind stencil point!');
        end
        if FM==0
            if single(e1+fc22(1,4)+fc22(1,1))~=single(e1+dis(fc17(:,1),fc17(:,4)))
                error('Wrong calculation of futher down/upwind stencil point distance!');
            end
            if single(e1+fc22(1,2)+fc22(1,3))~=single(e1+dis(fc17(:,2),fc17(:,3)))
                error('Wrong calculation of down/upwind stencil point distance!');
            end
        else
            if single(e+fc22(1,4)+fc22(1,1))~=single(e+dis(fc17(:,1),fc17(:,4)))
                error('Wrong calculation of futher down/upwind stencil point distance!');
            end
            if single(e+fc22(1,2)+fc22(1,3))~=single(e+dis(fc17(:,2),fc17(:,3)))
                error('Wrong calculation of down/upwind stencil point distance!');
            end
        end
    end
    
    %% FC{23}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc18=FC{18};
        fc23=zeros(1,length(fc16));
        for i=1:length(fc16)
            if fc16(1,i)==0
                if sum(fc18(:,i))==0
                    error('Boundary should be acrossed by the stencil!');
                end
                Nd1=NODE{fc18(1,i)};
                Nd2=NODE{fc18(2,i)};
                C_mid=(Nd1{3}+Nd2{3})/2;
                on_top=on_edge(C_mid,TL,TR);
                on_right=on_edge(C_mid,TR,BR);
                on_bottom=on_edge(C_mid,BL,BR);
                on_left=on_edge(C_mid,BL,TL);
                if (on_top+on_right+on_bottom+on_left)~=0 && (on_top+on_right+on_bottom+on_left)~=1
                    error('The boundary intercpt points are not on boundary!');
                end
                if on_top
                    fc23(1,i)=1;
                elseif on_right
                    fc23(1,i)=2;
                elseif on_bottom
                    fc23(1,i)=3;
                elseif on_left
                    fc23(1,i)=4;
                else
                    fc23(1,i)=-1; % On inner boundary
                end
            end
        end
        FC{23}=fc23;
        FACE{l}=FC;
    end
    % check FC{23}
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc17=FC{17};
        fc23=FC{23};
        for i=1:length(fc16)
            if fc16(1,i)==0
                if fc23(1,i)==0
                    if (fc17(1,i)<X1 || fc17(1,i)>X2) || (fc17(2,i)<Y1 || fc17(2,i)>Y2)
                        error('The stencil point should be within the outer boundaries!');
                    end
                elseif fc23(1,i)==1
%                     if (fc17(1,i)<X1 || fc17(1,i)>X2) || (fc17(2,i)<Y2)
                    if (fc17(2,i)<Y2)
                        error('The stencil point should be out of the top boundaries!');
                    end
                elseif fc23(1,i)==2
%                     if (fc17(1,i)<X2) || (fc17(2,i)<Y1 || fc17(2,i)>Y2)
                    if (fc17(1,i)<X2)
                        error('The stencil point should be out of the right boundaries!');
                    end
                elseif fc23(1,i)==3
%                     if (fc17(1,i)<X1 || fc17(1,i)>X2) || (fc17(2,i)>Y1)
                    if (fc17(2,i)>Y1)
                        error('The stencil point should be out of the bottom boundaries!');
                    end
                elseif fc23(1,i)==4
%                     if (fc17(1,i)>X1) || (fc17(2,i)<Y1 || fc17(2,i)>Y2)
                    if (fc17(1,i)>X1)
                        error('The stencil point should be out of the left boundaries!');
                    end
                elseif fc23(1,i)==-1
                    if (fc17(1,i)<X1 || fc17(1,i)>X2) || (fc17(2,i)<Y1 || fc17(2,i)>Y2)
                        error('The stencil point should be within the outer boundaries!');
                    end
                else
                    error('Wrong flag for boundary acrossed!');
                end
            else
                if (fc17(1,i)<X1 || fc17(1,i)>X2) || (fc17(2,i)<Y1 || fc17(2,i)>Y2)
                    error('The stencil point should be within the outer boundaries!');
                end
            end
        end
    end
    
    %% FC{24}
    for l=1:O
        FC=FACE{l};
        fc24=0;
        if FC{2}>0 % Type 2
            fc24=2;
        elseif FC{2}<0 % Type 1
            fc24=1;
        else % Type 1,3 and 4
            fc16=FC{16};
            fc20=FC{20};
            fc23=FC{23};
            if length(union(0,fc16))<=length(fc16) && (fc23(1,1)>=0 && fc23(1,end)>=0)% Face type 3
                fc24=3;
            else % Face type 1 and 4
                fc24=1;
                for i=1:length(fc16)
                    if fc20(1,i)==1 || fc20(1,i)==2
                        if fc20(1,i)==2
                            Nd1=NODE{fc20(2,i)};
                            Nd2=NODE{fc20(3,i)};
                            if Nd1{2}>0 || Nd2{2}>0
                                fc24=4;
                                break;
                            end
                        else
                            Nd=NODE{fc20(2,i)};
                            if Nd{2}>0
                                fc24=4;
                                break;
                            end
                        end
                    end
                end
            end
        end
        FC{24}=fc24;
        FACE{l}=FC;
    end
    % Check FC{24}
    counter1=0;
    counter2=0;
    counter3=0;
    for l=1:O
        FC=FACE{l};
        fc24=FC{24};
        fc16=FC{16};
        fc23=FC{23};
        if fc24==1
            if FC{2}<0
                counter1=counter1+1;
            elseif FC{2}==0
                counter3=counter3+1;
                if length(union(0,fc16))>length(fc16)
                    if sum(fc23)~=0
                        error('Logic error!');
                    end
                elseif length(union(0,fc16))==length(fc16)
                    if fc23(1,1)>=0 && fc23(1,end)>=0
                        error('The stencil of Type 1 face could only intercept inner boundary!');
                    end
                else
                    error('Logic error!')
                end
            else
                error('Type 1 face could be outer boundary!')
            end
        elseif fc24==2
            if FC{2}<=0
                error('Wrong value for face type 2!');
            end
            counter2=counter2+1;
        elseif fc24==3
            counter3=counter3+1;
            if FC{2}~=0
                error('Wrong value for face type 3!');
            end
            if length(union(0,fc16))==length(fc16)
                if fc23(1,1)<=0 && fc23(1,end)<=0
                    error('The stencil of Type 3 face could only intercept outer boundary!');
                end
            elseif length(union(0,fc16))<length(fc16)
                if fc23(1,1)<=0 || fc23(1,end)<=0
                    error('The stencil of Type 3 face could only intercept outer boundary!');
                end
            else
                error('Logic error!')
            end
        elseif fc24==4
            if FC{2}~=0
                error('Wrong value for face type 3!');
            end
            counter3=counter3+1;
        else
            error('Wrong value for FC{24}!');
        end
    end
    if counter1+counter2+counter3~=O
        error('Logic error!');
    end
    % Check whether type-3 face crosses the inner boundary
    for l=1:O
        FC=FACE{l};
        if FC{24}==3
            fc23=FC{23};
            for i=1:length(fc23)
                if fc23(1,i)<0
                    error('Refine the mesh between the inner and outer boundaries');
                end
            end
        end
    end
    % Check whether type-4 face crosses the inner boundary
    for l=1:O
        FC=FACE{l};
        if FC{24}==4
            fc23=FC{23};
            for i=1:length(fc23)
                if fc23(1,i)~=0
                    error('Refine the mesh between the inner and outer boundaries');
                end
            end
        end
    end
    % Check whether type-1 face voilates its defination
    for l=1:O
        FC=FACE{l};
        if FC{24}==1
            fc23=FC{23};
            fc20=FC{20};
            for i=1:length(fc23)
                if fc23(1,i)>0
                    error('All stencil points of type-1 face should be within in the outer boundaries!');
                else
                    if fc20(1,i)==1 || fc20(1,i)==2
                        if fc20(1,i)==2
                            Nd1=NODE{fc20(2,i)};
                            Nd2=NODE{fc20(3,i)};
                            if Nd1{2}>0 || Nd2{2}>0
                                error('All enclosing points for the stencil points of type-1 face should be either cell centroids or inner boundary nodes!');
                            end
                        else
                            Nd=NODE{fc20(2,i)};
                            if Nd{2}>0
                                error('All enclosing points for the stencil points of type-1 face should be either cell centroids or inner boundary nodes!');
                            end
                        end
                    end
                end
            end
        end
    end
end


%% FC{25~30}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || (FC{24}==3 || FC{24}==4)
        break;
    end
end
if (isempty(FC{25}) || isempty(FC{26})) || (isempty(FC{27}) || isempty(FC{28})) || (isempty(FC{29}) || isempty(FC{30}))
    %% FC{25}, FC{26} and FC{27}
    counter1=0;
    counter3=0;
    for l=1:O
        FC=FACE{l};
        fc24=FC{24};
        if fc24==2 || fc24==3 % Face type 2 and 3
            fc16=FC{16};
            fc23=FC{23};
            C_b=FC{7};
            counter1=counter1+1;
            fc17=FC{17};
            fc18=FC{18};
            fc22=FC{22};
            fc25=fc16;
            fc26=fc17;
            fc27=fc22;
            if ((fc16(1,1)==0 && fc16(1,2)==0) && (fc16(1,3)~=0 && fc16(1,4)~=0)) || ((fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)==0 && fc16(1,4)==0)) % Outer boundary face
                if (fc16(1,1)==0 && fc16(1,2)==0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    if fc23(1,1)~=fc23(1,2)
                        error('Error in FC{23}!');
                    end
                    if sum(fc18(:,1)~=fc18(:,2))>0
                        error('Error in FC{18}');
                    end
                    Nd1=NODE{fc18(1,1)};
                    Nd2=NODE{fc18(2,1)};
                    L=dis(Nd1{3},Nd2{3});
                    v=Nd2{3}-Nd1{3};
                    v=v/norm(v);
                    C_i=Nd1{3}+fc18(3,1)*L*v;
                    if fc23(1,1)==1
                        C_i_m=C_i+[0;-(Y2-Y1)];
                        if ~on_edge(C_i_m,BL,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==2
                        C_i_m=C_i+[-(X2-X1);0];
                        if ~on_edge(C_i_m,BL,TL);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==3
                        C_i_m=C_i+[0;(Y2-Y1)];
                        if ~on_edge(C_i_m,TL,TR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==4
                        C_i_m=C_i+[(X2-X1);0];
                        if ~on_edge(C_i_m,TR,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    else
                        error('Logic error!');
                    end
                    neigh_up_m=FC{14};
                    neigh_down_m=FC{15};
                    cell_m=setxor(fc16(1,3),[neigh_up_m(1,1),neigh_down_m(1,1)]);
                    if length(cell_m)~=1
                        error('The mirror cell is not found!');
                    end
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)==0 && fc16(1,4)==0)
                    if fc23(1,3)~=fc23(1,4)
                        error('Logic error!');
                    end
                    if sum(fc18(:,3)~=fc18(:,4))>0
                        error('Error in FC{18}');
                    end
                    Nd1=NODE{fc18(1,4)};
                    Nd2=NODE{fc18(2,4)};
                    L=dis(Nd1{3},Nd2{3});
                    v=Nd2{3}-Nd1{3};
                    v=v/norm(v);
                    C_i=Nd1{3}+fc18(3,4)*L*v;
                    if fc23(1,4)==1
                        C_i_m=C_i+[0;-(Y2-Y1)];
                        if ~on_edge(C_i_m,BL,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==2
                        C_i_m=C_i+[-(X2-X1);0];
                        if ~on_edge(C_i_m,BL,TL);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==3
                        C_i_m=C_i+[0;(Y2-Y1)];
                        if ~on_edge(C_i_m,TL,TR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==4
                        C_i_m=C_i+[(X2-X1);0];
                        if ~on_edge(C_i_m,TR,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    else
                        error('Logic error!');
                    end
                    neigh_up_m=FC{14};
                    neigh_down_m=FC{15};
                    cell_m=setxor(fc16(1,2),[neigh_up_m(1,1),neigh_down_m(1,1)]);
                    if length(cell_m)~=1
                        error('The mirror cell is not found!');
                    end
                else
                    error('Logic error!');
                end
                if single(e+dis(C_i,C_b))>single(e)
                    error('Wrong intercept point!');
                end
                %% Find the coordinate of close stencil point
                Cell_m=CELL{cell_m};
                Nd1=NODE{Cell_m{7}};
                Nd2=NODE{Cell_m{8}};
                Nd3=NODE{Cell_m{9}};
                C_j_c=norm_joint(C_i_m,FC{4}',Cell_m{5});
                if ~in_triangle(C_j_c,Nd1{3},Nd2{3},Nd3{3})
                    error('The generated point is not within the current cell!');
                end
                %% Find the coordinate of further stencil point
                Nd1=NODE{Cell_m{7}};
                Nd2=NODE{Cell_m{8}};
                Nd3=NODE{Cell_m{9}};
                if on_edge(C_i_m,Nd1{3},Nd2{3})
                    Nd_opp_face_f=Nd3;
                elseif on_edge(C_i_m,Nd2{3},Nd3{3})
                    Nd_opp_face_f=Nd1;
                elseif on_edge(C_i_m,Nd3{3},Nd1{3})
                    Nd_opp_face_f=Nd2;
                else
                    error('Logic error!');
                end
                
                cell_group_f=setxor(Cell_m{1},Nd_opp_face_f{5});
                if (Nd_opp_face_f{4}-length(cell_group_f))~=1
                    error('The upwind cell should be excluded!');
                end
                c_j_f=zeros(2,Nd_opp_face_f{4}-1);
                f_found_counter=0;
                f_found_cell=0;
                c_j_f_found=zeros(2,1);
                for i=1:Nd_opp_face_f{4}-1
                    P=CELL{cell_group_f(i)};
                    Nd1=NODE{P{7}};
                    Nd2=NODE{P{8}};
                    Nd3=NODE{P{9}};
                    c_j_f(:,i)=norm_joint(C_i_m,FC{4}',P{5});
                    if in_triangle(c_j_f(:,i),Nd1{3},Nd2{3},Nd3{3})
                        f_found_counter=f_found_counter+1;
                        f_found_cell(f_found_counter)=cell_group_f(i);
                        c_j_f_found(:,f_found_counter)=c_j_f(:,i);
                    end
                end
                if f_found_counter==0
                    error('The further cell is not found!');
                elseif f_found_counter==1
                    if f_found_cell==0
                        error('Logic error!');
                    end
                    cell_j_f=f_found_cell;
                    C_j_f=c_j_f_found;
                else % More than one further upwind cell is found
                    Dist=zeros(1,f_found_counter);
                    for i=1:f_found_counter
                        Dist(i)=dis(C_i_m,c_j_f_found(:,i));
                    end
                    Dist_min=min(Dist);
                    for i=1:f_found_counter
                        if single(e+Dist(i))==single(e+Dist_min)
                            break;
                        end
                    end
                    if i==f_found_counter
                        if single(e+Dist(i))~=single(e+Dist_min)
                            error('Logic error!');
                        end
                    end
                    cell_j_f=f_found_cell(i);
                    C_j_f=c_j_f_found(:,i);
                end
                % Check
                cell_found=CELL{cell_j_f};
                Nd1=NODE{cell_found{7}};
                Nd2=NODE{cell_found{8}};
                Nd3=NODE{cell_found{9}};
                if ~in_triangle(C_j_f,Nd1{3},Nd2{3},Nd3{3})
                    error('The further stencil point is not within the found cell!');
                end
                if dis(C_i_m,C_j_f)<dis(C_i_m,C_j_c)
                    error('The further stencil should be further away from the face!');
                end
                %% fill
                if (fc16(1,1)==0 && fc16(1,2)==0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    fc25(1,1)=cell_j_f;
                    fc26(:,1)=C_j_f;
                    fc27(1,1)=dis(C_i_m,C_j_f);
                    
                    fc25(1,2)=cell_m;
                    fc26(:,2)=C_j_c;
                    fc27(1,2)=dis(C_i_m,C_j_c);
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)==0 && fc16(1,4)==0)
                    fc25(1,4)=cell_j_f;
                    fc26(:,4)=C_j_f;
                    fc27(1,4)=dis(C_i_m,C_j_f);
                    
                    fc25(1,3)=cell_m;
                    fc26(:,3)=C_j_c;
                    fc27(1,3)=dis(C_i_m,C_j_c);
                else
                    error('Logic error!');
                end
                %% Final check
                if (fc16(1,1)==0 && fc16(1,2)==0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    c_dd=fc26(:,1);
                    c_d=fc26(:,2);
                    if fc23(1,1)==1
                        c_dd=c_dd+[0;(Y2-Y1)];
                        c_d=c_d+[0;(Y2-Y1)];
                    elseif fc23(1,1)==2
                        c_dd=c_dd+[(X2-X1);0];
                        c_d=c_d+[(X2-X1);0];
                    elseif fc23(1,1)==3
                        c_dd=c_dd+[0;-(Y2-Y1)];
                        c_d=c_d+[0;-(Y2-Y1)];
                    elseif fc23(1,1)==4
                        c_dd=c_dd+[-(X2-X1);0];
                        c_d=c_d+[-(X2-X1);0];
                    else
                        error('Logic error!');
                    end
                    if FM==0
                        if (single(e1+dis(c_dd,C_b))~=single(e1+fc27(1,1)) || single(e1+dis(c_d,C_b))~=single(e1+fc27(1,2))) || (single(e1+dis(c_dd,fc26(:,4)))~=single(e1+fc27(1,1)+fc27(1,4)) || single(e1+dis(c_d,fc26(:,3)))~=single(e1+fc27(1,2)+fc27(1,3)))
                            error('FC{27} contains false info!');
                        end
                    else
                        if (single(e1+dis(c_dd,C_b))~=single(e1+fc27(1,1)) || single(e1+dis(c_d,C_b))~=single(e1+fc27(1,2))) || (single(e1+dis(c_dd,fc26(:,4)))~=single(e1+fc27(1,1)+fc27(1,4)) || single(e1+dis(c_d,fc26(:,3)))~=single(e1+fc27(1,2)+fc27(1,3)))
                            error('FC{27} contains false info!');
                        end
                    end
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)==0 && fc16(1,4)==0)
                    c_uu=fc26(:,4);
                    c_u=fc26(:,3);
                    if fc23(1,4)==1
                        c_uu=c_uu+[0;(Y2-Y1)];
                        c_u=c_u+[0;(Y2-Y1)];
                    elseif fc23(1,4)==2
                        c_uu=c_uu+[(X2-X1);0];
                        c_u=c_u+[(X2-X1);0];
                    elseif fc23(1,4)==3
                        c_uu=c_uu+[0;-(Y2-Y1)];
                        c_u=c_u+[0;-(Y2-Y1)];
                    elseif fc23(1,4)==4
                        c_uu=c_uu+[-(X2-X1);0];
                        c_u=c_u+[-(X2-X1);0];
                    else
                        error('Logic error!');
                    end
                    if FM==0
                        if (single(e1+dis(c_uu,C_b))~=single(e1+fc27(1,4)) || single(e1+dis(c_u,C_b))~=single(e1+fc27(1,3))) || (single(e1+dis(c_uu,fc26(:,1)))~=single(e1+fc27(1,1)+fc27(1,4)) || single(e1+dis(c_u,fc26(:,2)))~=single(e1+fc27(1,2)+fc27(1,3)))
                            error('FC{27} contains false info!');
                        end
                    else
                        if (single(e1+dis(c_uu,C_b))~=single(e1+fc27(1,4)) || single(e1+dis(c_u,C_b))~=single(e1+fc27(1,3))) || (single(e1+dis(c_uu,fc26(:,1)))~=single(e1+fc27(1,1)+fc27(1,4)) || single(e1+dis(c_u,fc26(:,2)))~=single(e1+fc27(1,2)+fc27(1,3)))
                            error('FC{27} contains false info!');
                        end
                    end
                else
                    error('Logic error!');
                end
            elseif ((fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0)) || ((fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0)) % Interior face that has either further downwind or further upwind stencil point out of outer boundary
                if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    Nd1=NODE{fc18(1,1)};
                    Nd2=NODE{fc18(2,1)};
                    L=dis(Nd1{3},Nd2{3});
                    v=Nd2{3}-Nd1{3};
                    v=v/norm(v);
                    C_i=Nd1{3}+fc18(3,1)*L*v;
                    if fc23(1,1)==1
                        C_i_m=C_i+[0;-(Y2-Y1)];
                        if ~on_edge(C_i_m,BL,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==2
                        C_i_m=C_i+[-(X2-X1);0];
                        if ~on_edge(C_i_m,BL,TL);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==3
                        C_i_m=C_i+[0;(Y2-Y1)];
                        if ~on_edge(C_i_m,TL,TR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,1)==4
                        C_i_m=C_i+[(X2-X1);0];
                        if ~on_edge(C_i_m,TR,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    else
                        error('Logic error!');
                    end
                    FCC=intersect(Nd1{13},Nd2{13});
                    if length(FCC)~=1
                        error('Common face is not found!');
                    end
                    FCC=FACE{FCC};
                    if FCC{2}<=0
                        error('The common face should be outer boundary face!');
                    end
                    neigh_up_nm=FCC{12};
                    neigh_down_nm=FCC{13};
                    neigh_up_m=FCC{14};
                    neigh_down_m=FCC{15};
                    cell_m=setxor(setxor(0,[neigh_up_nm(1,1),neigh_down_nm(1,1)]),[neigh_up_m(1,1),neigh_down_m(1,1)]);
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0)
                    Nd1=NODE{fc18(1,4)};
                    Nd2=NODE{fc18(2,4)};
                    L=dis(Nd1{3},Nd2{3});
                    v=Nd2{3}-Nd1{3};
                    v=v/norm(v);
                    C_i=Nd1{3}+fc18(3,4)*L*v;
                    if fc23(1,4)==1
                        C_i_m=C_i+[0;-(Y2-Y1)];
                        if ~on_edge(C_i_m,BL,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==2
                        C_i_m=C_i+[-(X2-X1);0];
                        if ~on_edge(C_i_m,BL,TL);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==3
                        C_i_m=C_i+[0;(Y2-Y1)];
                        if ~on_edge(C_i_m,TL,TR);
                            error('The mirror intercept point is wrong!');
                        end
                    elseif fc23(1,4)==4
                        C_i_m=C_i+[(X2-X1);0];
                        if ~on_edge(C_i_m,TR,BR);
                            error('The mirror intercept point is wrong!');
                        end
                    else
                        error('Logic error!');
                    end
                    FCC=intersect(Nd1{13},Nd2{13});
                    if length(FCC)~=1
                        error('Common face is not found!');
                    end
                    FCC=FACE{FCC};
                    if FCC{2}<=0
                        error('The common face should be outer boundary face!');
                    end
                    neigh_up_nm=FCC{12};
                    neigh_down_nm=FCC{13};
                    neigh_up_m=FCC{14};
                    neigh_down_m=FCC{15};
                    cell_m=setxor(setxor(0,[neigh_up_nm(1,1),neigh_down_nm(1,1)]),[neigh_up_m(1,1),neigh_down_m(1,1)]);
                    if length(cell_m)~=1
                        error('The mirror cell is not found!');
                    end
                else
                    error('Logic error!');
                end
                
                %%
                Cell_m=CELL{cell_m};
                Nd1=NODE{Cell_m{7}};
                Nd2=NODE{Cell_m{8}};
                Nd3=NODE{Cell_m{9}};
                
                cell_group_f=union(Nd1{5},union(Nd2{5},Nd3{5}));
                c_j_f=zeros(2,length(cell_group_f));
                f_found_counter=0;
                f_found_cell=0;
                c_j_f_found=zeros(2,1);
                for i=1:length(cell_group_f)
                    P=CELL{cell_group_f(i)};
                    Nd1=NODE{P{7}};
                    Nd2=NODE{P{8}};
                    Nd3=NODE{P{9}};
                    c_j_f(:,i)=norm_joint(C_i_m,FC{4}',P{5});
                    if in_triangle(c_j_f(:,i),Nd1{3},Nd2{3},Nd3{3})
                        f_found_counter=f_found_counter+1;
                        f_found_cell(f_found_counter)=cell_group_f(i);
                        c_j_f_found(:,f_found_counter)=c_j_f(:,i);
                    end
                end
                if f_found_counter==0
                    error('The further cell is not found!');
                elseif f_found_counter==1
                    if f_found_cell==0
                        error('Logic error!');
                    end
                    cell_j_f=f_found_cell;
                    C_j_f=c_j_f_found;
                else % More than one further upwind cell is found
                    Dist=zeros(1,f_found_counter);
                    for i=1:f_found_counter
                        Dist(i)=dis(C_i_m,c_j_f_found(:,i));
                    end
                    Dist_min=min(Dist);
                    for i=1:f_found_counter
                        if single(e+Dist(i))==single(e+Dist_min)
                            break;
                        end
                    end
                    if i==f_found_counter
                        if single(e+Dist(i))~=single(e+Dist_min)
                            error('Logic error!');
                        end
                    end
                    cell_j_f=f_found_cell(i);
                    C_j_f=c_j_f_found(:,i);
                end
                %% fill
                if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    fc25(1,1)=cell_j_f;
                    fc26(:,1)=C_j_f;
                    fc27(1,1)=dis(C_i_m,C_j_f)+fc22(1,2)+dis(C_i,fc17(:,2));
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0)
                    fc25(1,4)=cell_j_f;
                    fc26(:,4)=C_j_f;
                    fc27(1,4)=dis(C_i_m,C_j_f)+fc22(1,3)+dis(C_i,fc17(:,3));
                else
                    error('Logic error!');
                end
                %% Final check
                if (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)~=0)
                    c_dd=fc26(:,1);
                    if fc23(1,1)==1
                        c_dd=c_dd+[0;(Y2-Y1)];
                    elseif fc23(1,1)==2
                        c_dd=c_dd+[(X2-X1);0];
                    elseif fc23(1,1)==3
                        c_dd=c_dd+[0;-(Y2-Y1)];
                    elseif fc23(1,1)==4
                        c_dd=c_dd+[-(X2-X1);0];
                    else
                        error('Logic error!');
                    end
                    if FM==0
                        if single(e1+dis(c_dd,C_b))~=single(e1+fc27(1,1)) || single(e1+dis(c_dd,fc26(:,4)))~=single(e1+fc27(1,1)+fc27(1,4))
                            error('FC{27} contains false info!');
                        end
                    else
                        if single(e+dis(c_dd,C_b))~=single(e+fc27(1,1)) || single(e+dis(c_dd,fc26(:,4)))~=single(e+fc27(1,1)+fc27(1,4))
                            error('FC{27} contains false info!');
                        end
                    end
                elseif (fc16(1,1)~=0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0)
                    c_uu=fc26(:,4);
                    if fc23(1,4)==1
                        c_uu=c_uu+[0;(Y2-Y1)];
                    elseif fc23(1,4)==2
                        c_uu=c_uu+[(X2-X1);0];
                    elseif fc23(1,4)==3
                        c_uu=c_uu+[0;-(Y2-Y1)];
                    elseif fc23(1,4)==4
                        c_uu=c_uu+[-(X2-X1);0];
                    else
                        error('Logic error!');
                    end
                    if FM==0
                        if single(e1+dis(c_uu,C_b))~=single(e1+fc27(1,4))  || single(e1+dis(c_uu,fc26(:,1)))~=single(e1+fc27(1,1)+fc27(1,4))
                            error('FC{27} contains false info!');
                        end
                    else
                        if single(e+dis(c_uu,C_b))~=single(e+fc27(1,4))  || single(e+dis(c_uu,fc26(:,1)))~=single(e+fc27(1,1)+fc27(1,4))
                            error('FC{27} contains false info!');
                        end
                    end
                else
                    error('Logic error!');
                end
            elseif (fc16(1,1)==0 && fc16(1,2)~=0) && (fc16(1,3)~=0 && fc16(1,4)==0) % Interior face that both further downwind and further upwind stencil points out of outer boundary
                %% First intercept point
                Nd1=NODE{fc18(1,1)};
                Nd2=NODE{fc18(2,1)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i_1=Nd1{3}+fc18(3,1)*L*v;
                if fc23(1,1)==1
                    C_i_m_1=C_i_1+[0;-(Y2-Y1)];
                    if ~on_edge(C_i_m_1,BL,BR);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,1)==2
                    C_i_m_1=C_i_1+[-(X2-X1);0];
                    if ~on_edge(C_i_m_1,BL,TL);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,1)==3
                    C_i_m_1=C_i_1+[0;(Y2-Y1)];
                    if ~on_edge(C_i_m_1,TL,TR);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,1)==4
                    C_i_m_1=C_i_1+[(X2-X1);0];
                    if ~on_edge(C_i_m_1,TR,BR);
                        error('The mirror intercept point is wrong!');
                    end
                else
                    error('Logic error!');
                end
                FCC=intersect(Nd1{13},Nd2{13});
                if length(FCC)~=1
                    error('Common face is not found!');
                end
                FCC=FACE{FCC};
                if FCC{2}<=0
                    error('The common face should be outer boundary face!');
                end
                neigh_up_nm_1=FCC{12};
                neigh_down_nm_1=FCC{13};
                neigh_up_m_1=FCC{14};
                neigh_down_m_1=FCC{15};
                cell_m_1=setxor(setxor(0,[neigh_up_nm_1(1,1),neigh_down_nm_1(1,1)]),[neigh_up_m_1(1,1),neigh_down_m_1(1,1)]);
                if length(cell_m_1)~=1
                    error('The mirror cell is not found!');
                end
                Cell_m_1=CELL{cell_m_1};
                
                Nd1=NODE{Cell_m_1{7}};
                Nd2=NODE{Cell_m_1{8}};
                Nd3=NODE{Cell_m_1{9}};
                
                cell_group_f=union(Nd1{5},union(Nd2{5},Nd3{5}));
                c_j_f=zeros(2,length(cell_group_f));
                f_found_counter=0;
                f_found_cell=0;
                c_j_f_found=zeros(2,1);
                for i=1:length(cell_group_f)
                    P=CELL{cell_group_f(i)};
                    Nd1=NODE{P{7}};
                    Nd2=NODE{P{8}};
                    Nd3=NODE{P{9}};
                    c_j_f(:,i)=norm_joint(C_i_m_1,FC{4}',P{5});
                    if in_triangle(c_j_f(:,i),Nd1{3},Nd2{3},Nd3{3})
                        f_found_counter=f_found_counter+1;
                        f_found_cell(f_found_counter)=cell_group_f(i);
                        c_j_f_found(:,f_found_counter)=c_j_f(:,i);
                    end
                end
                if f_found_counter==0
                    error('The further cell is not found!');
                elseif f_found_counter==1
                    if f_found_cell==0
                        error('Logic error!');
                    end
                    cell_j_f_1=f_found_cell;
                    C_j_f_1=c_j_f_found;
                else % More than one further upwind cell is found
                    Dist=zeros(1,f_found_counter);
                    for i=1:f_found_counter
                        Dist(i)=dis(C_i_m_1,c_j_f_found(:,i));
                    end
                    Dist_min=min(Dist);
                    for i=1:f_found_counter
                        if single(e+Dist(i))==single(e+Dist_min)
                            break;
                        end
                    end
                    if i==f_found_counter
                        if single(e+Dist(i))~=single(e+Dist_min)
                            error('Logic error!');
                        end
                    end
                    cell_j_f_1=f_found_cell(i);
                    C_j_f_1=c_j_f_found(:,i);
                end
                % Check
                Cell_m_1=CELL{cell_j_f_1};
                Nd1=NODE{Cell_m_1{7}};
                Nd2=NODE{Cell_m_1{8}};
                Nd3=NODE{Cell_m_1{9}};
                if ~in_triangle(C_j_f_1,Nd1{3},Nd2{3},Nd3{3})
                    error('The generated point is not within the current cell!');
                end
                
                
                %% Second intercept point
                Nd1=NODE{fc18(1,4)};
                Nd2=NODE{fc18(2,4)};
                L=dis(Nd1{3},Nd2{3});
                v=Nd2{3}-Nd1{3};
                v=v/norm(v);
                C_i_2=Nd1{3}+fc18(3,4)*L*v;
                if fc23(1,4)==1
                    C_i_m_2=C_i_2+[0;-(Y2-Y1)];
                    if ~on_edge(C_i_m_2,BL,BR);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,4)==2
                    C_i_m_2=C_i_2+[-(X2-X1);0];
                    if ~on_edge(C_i_m_2,BL,TL);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,4)==3
                    C_i_m_2=C_i_2+[0;(Y2-Y1)];
                    if ~on_edge(C_i_m_2,TL,TR);
                        error('The mirror intercept point is wrong!');
                    end
                elseif fc23(1,4)==4
                    C_i_m_2=C_i_2+[(X2-X1);0];
                    if ~on_edge(C_i_m_2,TR,BR);
                        error('The mirror intercept point is wrong!');
                    end
                else
                    error('Logic error!');
                end
                FCC=intersect(Nd1{13},Nd2{13});
                if length(FCC)~=1
                    error('Common face is not found!');
                end
                FCC=FACE{FCC};
                if FCC{2}<=0
                    error('The common face should be outer boundary face!');
                end
                neigh_up_nm_2=FCC{12};
                neigh_down_nm_2=FCC{13};
                neigh_up_m_2=FCC{14};
                neigh_down_m_2=FCC{15};
                cell_m_2=setxor(setxor(0,[neigh_up_nm_2(1,1),neigh_down_nm_2(1,1)]),[neigh_up_m_2(1,1),neigh_down_m_2(1,1)]);
                if length(cell_m_2)~=1
                    error('The mirror cell is not found!');
                end
                Cell_m_2=CELL{cell_m_2};
                
                Nd1=NODE{Cell_m_2{7}};
                Nd2=NODE{Cell_m_2{8}};
                Nd3=NODE{Cell_m_2{9}};
                
                cell_group_f=union(Nd1{5},union(Nd2{5},Nd3{5}));
                c_j_f=zeros(2,length(cell_group_f));
                f_found_counter=0;
                f_found_cell=0;
                c_j_f_found=zeros(2,1);
                for i=1:length(cell_group_f)
                    P=CELL{cell_group_f(i)};
                    Nd1=NODE{P{7}};
                    Nd2=NODE{P{8}};
                    Nd3=NODE{P{9}};
                    c_j_f(:,i)=norm_joint(C_i_m_2,FC{4}',P{5});
                    if in_triangle(c_j_f(:,i),Nd1{3},Nd2{3},Nd3{3})
                        f_found_counter=f_found_counter+1;
                        f_found_cell(f_found_counter)=cell_group_f(i);
                        c_j_f_found(:,f_found_counter)=c_j_f(:,i);
                    end
                end
                if f_found_counter==0
                    error('The further cell is not found!');
                elseif f_found_counter==1
                    if f_found_cell==0
                        error('Logic error!');
                    end
                    cell_j_f_2=f_found_cell;
                    C_j_f_2=c_j_f_found;
                else % More than one further upwind cell is found
                    Dist=zeros(1,f_found_counter);
                    for i=1:f_found_counter
                        Dist(i)=dis(C_i_m_2,c_j_f_found(:,i));
                    end
                    Dist_min=min(Dist);
                    for i=1:f_found_counter
                        if single(e+Dist(i))==single(e+Dist_min)
                            break;
                        end
                    end
                    if i==f_found_counter
                        if single(e+Dist(i))~=single(e+Dist_min)
                            error('Logic error!');
                        end
                    end
                    cell_j_f_2=f_found_cell(i);
                    C_j_f_2=c_j_f_found(:,i);
                end
                % Check
                Cell_m_2=CELL{cell_j_f_2};
                Nd1=NODE{Cell_m_2{7}};
                Nd2=NODE{Cell_m_2{8}};
                Nd3=NODE{Cell_m_2{9}};
                if ~in_triangle(C_j_f_2,Nd1{3},Nd2{3},Nd3{3})
                    error('The generated point is not within the current cell!');
                end
                
                %% fill
                fc25(1,1)=cell_j_f_1;
                fc26(:,1)=C_j_f_1;
                fc27(1,1)=dis(C_i_m_1,C_j_f_1)+fc22(1,2)+dis(C_i_1,fc17(:,2));
                
                fc25(1,4)=cell_j_f_2;
                fc26(:,4)=C_j_f_2;
                fc27(1,4)=dis(C_i_m_2,C_j_f_2)+fc22(1,3)+dis(C_i_2,fc17(:,3));
                
                %% Final check
                c_dd=fc26(:,1);
                if fc23(1,1)==1
                    c_dd=c_dd+[0;(Y2-Y1)];
                elseif fc23(1,1)==2
                    c_dd=c_dd+[(X2-X1);0];
                elseif fc23(1,1)==3
                    c_dd=c_dd+[0;-(Y2-Y1)];
                elseif fc23(1,1)==4
                    c_dd=c_dd+[-(X2-X1);0];
                else
                    error('Logic error!');
                end
                if single(e1+dis(c_dd,C_b))~=single(e1+fc27(1,1))
                    error('FC{27} contains false info!');
                end
                
                c_uu=fc26(:,4);
                if fc23(1,4)==1
                    c_uu=c_uu+[0;(Y2-Y1)];
                elseif fc23(1,4)==2
                    c_uu=c_uu+[(X2-X1);0];
                elseif fc23(1,4)==3
                    c_uu=c_uu+[0;-(Y2-Y1)];
                elseif fc23(1,4)==4
                    c_uu=c_uu+[-(X2-X1);0];
                else
                    error('Logic error!');
                end
                if FM==0
                    if single(e1+dis(c_uu,C_b))~=single(e1+fc27(1,4))
                        error('FC{27} contains false info!');
                    end
                    
                    if single(e1+dis(c_uu,c_dd))~=single(e1+fc27(1,1)+fc27(1,4))
                        error('FC{27} contains false info!');
                    end
                else
                    if single(e+dis(c_uu,C_b))~=single(e+fc27(1,4))
                        error('FC{27} contains false info!');
                    end
                    
                    if single(e+dis(c_uu,c_dd))~=single(e+fc27(1,1)+fc27(1,4))
                        error('FC{27} contains false info!');
                    end
                end
            else
                error('Logic error!');
            end
            FC{25}=fc25;
            FC{26}=fc26;
            FC{27}=fc27;
            FACE{l}=FC;
        elseif FC{24}==4
            counter3=counter3+1;
            FC{25}=FC{16};
            FC{26}=FC{17};
            FC{27}=FC{22};
            FACE{l}=FC;
        end
    end
    
    % Check totality
    counter2=0; % For face type 1
    for l=1:O
        FC=FACE{l};
        if FC{2}<0 % Exclude the inner boundary
            counter2=counter2+1;
        elseif FC{2}==0
            fc16=FC{16};
            fc23=FC{23};
            if length(union(0,fc16))>length(fc16)
                fc20=FC{20};
                on_boundary_enclosing_node=0;
                for i=1:length(fc16)
                    if fc20(1,i)==1 || fc20(1,i)==2
                        if fc20(1,i)==2
                            Nd1=NODE{fc20(2,i)};
                            Nd2=NODE{fc20(3,i)};
                            if Nd1{2}>0 || Nd2{2}>0
                                on_boundary_enclosing_node=1;
                                break;
                            end
                        else
                            Nd=NODE{fc20(2,i)};
                            if Nd{2}>0
                                on_boundary_enclosing_node=1;
                                break;
                            end
                        end
                    end
                end
                if ~on_boundary_enclosing_node
                    counter2=counter2+1;
                end
            else
                if fc23(1,1)<0 || fc23(1,end)<0
                    counter2=counter2+1;
                end
            end
        end
    end
    if counter1+counter2+counter3~=O
        error('Unity is not formed!')
    end
    % Check FC{25}, FC{26}, FC{27}
    for l=1:O
        FC=FACE{l};
        c_b=FC{7};
        fc23=FC{23};
        fc24=FC{24};
        fc25=FC{25};
        fc26=FC{26};
        fc27=FC{27};
        % Check FC{25} and FC{26}
        for i=1:length(fc25)
            Cell_in=CELL{fc25(1,i)};
            Nd1=NODE{Cell_in{7}};
            Nd2=NODE{Cell_in{8}};
            Nd3=NODE{Cell_in{9}};
            if ~in_triangle(fc26(:,i),Nd1{3},Nd2{3},Nd3{3})
                error('FC{25} and/or FC{26} contain false info!');
            end
        end
        % Check FC{27}
        if fc24==2 || fc24==3
            for i=1:length(fc23)
                if fc23(1,i)==-1
                    ;
                elseif fc23(1,i)==0
                    if double(e1+fc27(1,i))~=double(e1+dis(c_b,fc26(:,i)))
                        error('FC{27} contains false info!');
                    end
                elseif fc23(1,i)==1
                    if single(e1+fc27(1,i))~=single(e1+dis(c_b,fc26(:,i)+[0;(Y2-Y1)]))
                        error('FC{27} contains false info!');
                    end
                elseif fc23(1,i)==2
                    if single(e1+fc27(1,i))~=single(e1+dis(c_b,fc26(:,i)+[(X2-X1);0]))
                        error('FC{27} contains false info!');
                    end
                elseif fc23(1,i)==3
                    if single(e1+fc27(1,i))~=single(e1+dis(c_b,fc26(:,i)-[0;(Y2-Y1)]))
                        error('FC{27} contains false info!');
                    end
                elseif fc23(1,i)==4
                    if single(e1+fc27(1,i))~=single(e1+dis(c_b,fc26(:,i)-[(X2-X1);0]))
                        error('FC{27} contains false info!');
                    end
                else
                    error('Logic error!');
                end
            end
        end
    end
    
    
    %% FC{28}
    e1=1e-7;
    for l=1:O
        FC=FACE{l};
        fc24=FC{24};
        if fc24==2 || fc24==3 % Face type 2 and 3
            fc16=FC{16};
            fc25=FC{25};
            fc26=FC{26};
            fc28=FC{19};
            for i=1:length(fc16)
                if fc16(1,i)==0
                    cell_in=CELL{fc25(1,i)};
                    Nd1=NODE{cell_in{7}};
                    Nd2=NODE{cell_in{8}};
                    Nd3=NODE{cell_in{9}};
                    C_c=cell_in{5};
                    if dis(fc26(:,i),C_c)<e1
                        in_zone_four=1;
                    else
                        in_zone_one=in_triangle(fc26(:,i),C_c,Nd1{3},Nd2{3});
                        in_zone_two=in_triangle(fc26(:,i),C_c,Nd2{3},Nd3{3});
                        in_zone_three=in_triangle(fc26(:,i),C_c,Nd3{3},Nd1{3});
                        in_zone_four=0;
                    end
                    % Check
                    if (in_zone_one+in_zone_two+in_zone_three)~=1 && (in_zone_one+in_zone_two+in_zone_three)~=0
                        error('The stencil point could no be in multiple zones!');
                    end
                    % Further determination
                    if in_zone_four
                        fc28(1,i)=4; % Located in zone 4
                    else
                        if in_zone_one
                            fc28(1,i)=1; % Located in zone 1
                        elseif in_zone_two
                            fc28(1,i)=2; % Located in zone 2
                        elseif in_zone_three
                            fc28(1,i)=3; % Located in zone 3
                        else
                            on_edge_one=on_edge(fc26(:,i),C_c,Nd1{3});
                            on_edge_two=on_edge(fc26(:,i),C_c,Nd2{3});
                            on_edge_three=on_edge(fc26(:,i),C_c,Nd3{3});
                            % Check
                            if (on_edge_one+on_edge_two+on_edge_three)~=1 && (on_edge_one+on_edge_two+on_edge_three)~=3
                                error('The stencil point could no be on multiple edges!');
                            end
                            if on_edge_one && (on_edge_two && on_edge_three)
                                % Check
                                if single(e+dis(C_c,fc26(:,i)))~=single(e)
                                    error('Logic error!');
                                end
                                fc28(1,i)=4;
                            elseif on_edge_one
                                fc28(1,i)=1;
                            elseif on_edge_two
                                fc28(1,i)=2;
                            elseif on_edge_three
                                fc28(1,i)=3;
                            else
                                error('The zone for the stencil point is not found!');
                            end
                        end
                    end
                end
            end
            FC{28}=fc28;
            %% Final data filling
            FACE{l}=FC;
        elseif FC{24}==4 % face type 4
            FC{28}=FC{19};
            %% Final data filling
            FACE{l}=FC;
        end
    end
    % Check for FC{28}
    for l=1:O
        FC=FACE{l};
        if FC{24}==2 || (FC{24}==3 || FC{24}==4)
            fc25=FC{25};
            fc26=FC{26};
            fc28=FC{28};
            for i=1:length(fc25)
                Cell_fix=CELL{fc25(1,i)};
                Nd1=NODE{Cell_fix{7}};
                Nd2=NODE{Cell_fix{8}};
                Nd3=NODE{Cell_fix{9}};
                C_fix=Cell_fix{5};
                if fc28(1,i)==1
                    C1=Nd1{3};
                    C2=Nd2{3};
                    if on_edge(fc26(:,i),C_fix,C1);
                        break;
                    end
                elseif fc28(1,i)==2
                    C1=Nd2{3};
                    C2=Nd3{3};
                    if on_edge(fc26(:,i),C_fix,C1);
                        break;
                    end
                elseif fc28(1,i)==3
                    C1=Nd3{3};
                    C2=Nd1{3};
                    if on_edge(fc26(:,i),C_fix,C1);
                        break;
                    end
                elseif fc28(1,i)==4
                    if dis(fc26(:,i),C_fix)>e1
                        error('The stencil point should be at the centroid!');
                    else
                        break;
                    end
                else
                    error('The points that enlose the stencil point cannot be all nodes!');
                end
                if ~in_triangle(fc26(:,i),C_fix,C1,C2)
                    error('Wrong zone ID!');
                end
            end
        end
    end
    
    %% FC{29}, {30}
    for l=1:O
        FC=FACE{l};
        fc24=FC{24};
        if fc24==2 || (fc24==3 || fc24==4) % Face type 2, 3 and 4
            fc25=FC{25};
            fc26=FC{26};
            fc28=FC{28};
            fc29=zeros(4,4);
            fc30=zeros(6,4);
            for i=1:length(fc25)
                if fc28(1,i)==4
                    fc29(:,i)=[0;fc25(1,i);fc25(1,i);fc25(1,i)];
                    fc30(:,i)=[fc26(:,i);fc26(:,i);fc26(:,i)];
                else
                    Cell_fix=CELL{fc25(1,i)};
                    C_fix=Cell_fix{5};
                    Face_target=FACE{Cell_fix{15+fc28(1,i)}};
                    Nd1=NODE{Face_target{8}};
                    Nd2=NODE{Face_target{9}};
                    if Face_target{10}~=0 && Face_target{11}~=0 % Both end nodes of target face are located on outer boundary
                        Onb1=on_which_boundary(Nd1{3},X1,X2,Y1,Y2);
                        Onb2=on_which_boundary(Nd2{3},X1,X2,Y1,Y2);
                        if (length(Onb1)==1 && length(Onb2)==2) || (length(Onb1)==2 && length(Onb2)==1)
                            if length(Onb1)==1 && length(Onb2)==2
                                Nd_b=Nd1;
                                Nd_c=Nd2;
                                Onb_b=Onb1;
                                Onb_c=Onb2;
                            elseif length(Onb1)==2 && length(Onb2)==1
                                Nd_b=Nd2;
                                Nd_c=Nd1;
                                Onb_b=Onb2;
                                Onb_c=Onb1;
                            else
                                error('Logic error!');
                            end
                            %% First centroid pool
                            cell_pool_1=setxor(Cell_fix{1},union(Nd_b{5},Nd_c{5})); % The centroid candidates on the local side of the periodic boundaries
                            Cell_pool_1=cell(1,length(cell_pool_1)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                            for s=1:length(cell_pool_1)
                                C=Cell_pool_1{s};
                                Cell_p=CELL{cell_pool_1(s)};
                                C{1,1}=Cell_p{1};
                                C{2,1}=Cell_p{5};
                                Cell_pool_1{s}=C;
                            end
                            %% Second centroid pool
                            % The sub pool from Nd_b
                            cell_pool_2_1=setxor(Nd_b{5},Nd_b{9}); % The centroid candidates on the periodic side of the periodic boundaries, where the boundary node is on regular boundary face
                            Cell_pool_2_1=cell(1,length(cell_pool_2_1)); % Temperary cell structure, which contains only the cell order number and centroid coordinates, where the boundary node is on regular boundary face
                            for s=1:length(cell_pool_2_1)
                                C=Cell_pool_2_1{s};
                                Cell_p=CELL{cell_pool_2_1(s)};
                                C{1,1}=Cell_p{1};
                                if Onb_b==1
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                elseif Onb_b==2
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                elseif Onb_b==3
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                elseif Onb_b==4
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                else
                                    error('Logic error!');
                                end
                                Cell_pool_2_1{s}=C;
                            end
                            % The sub pool from Nd_c
                            if length(union(Onb_c,[1,2]))==2 % Top right corner
                                Nd_c_1=NODE{N_I+N_L-1+N_H-2+N_L};
                                Nd_c_2=NODE{N_I+N_L-1+N_H-1};
                                Nd_c_3=NODE{N};
                                %
                                cell_pool_2_2_1=Nd_c_1{5};
                                Cell_pool_2_2_1=cell(1,length(cell_pool_2_2_1));
                                for s=1:length(cell_pool_2_2_1)
                                    C=Cell_pool_2_2_1{s};
                                    Cell_p=CELL{cell_pool_2_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);(Y2-Y1)];
                                    Cell_pool_2_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2_2=Nd_c_2{5};
                                Cell_pool_2_2_2=cell(1,length(cell_pool_2_2_2));
                                for s=1:length(cell_pool_2_2_2)
                                    C=Cell_pool_2_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                    Cell_pool_2_2_2{s}=C;
                                end
                                %
                                cell_pool_2_2_3=Nd_c_3{5};
                                Cell_pool_2_2_3=cell(1,length(cell_pool_2_2_3));
                                for s=1:length(cell_pool_2_2_3)
                                    C=Cell_pool_2_2_3{s};
                                    Cell_p=CELL{cell_pool_2_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                    Cell_pool_2_2_3{s}=C;
                                end
                            elseif length(union(Onb_c,[2,3]))==2 % Botom right corner
                                Nd_c_1=NODE{N};
                                Nd_c_2=NODE{N_I+N_L-1+N_H-2+N_L};
                                Nd_c_3=NODE{N_I+N_L-1};
                                %
                                cell_pool_2_2_1=Nd_c_1{5};
                                Cell_pool_2_2_1=cell(1,length(cell_pool_2_2_1));
                                for s=1:length(cell_pool_2_2_1)
                                    C=Cell_pool_2_2_1{s};
                                    Cell_p=CELL{cell_pool_2_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);-(Y2-Y1)];
                                    Cell_pool_2_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2_2=Nd_c_2{5};
                                Cell_pool_2_2_2=cell(1,length(cell_pool_2_2_2));
                                for s=1:length(cell_pool_2_2_2)
                                    C=Cell_pool_2_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                    Cell_pool_2_2_2{s}=C;
                                end
                                %
                                cell_pool_2_2_3=Nd_c_3{5};
                                Cell_pool_2_2_3=cell(1,length(cell_pool_2_2_3));
                                for s=1:length(cell_pool_2_2_3)
                                    C=Cell_pool_2_2_3{s};
                                    Cell_p=CELL{cell_pool_2_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                    Cell_pool_2_2_3{s}=C;
                                end
                            elseif length(union(Onb_c,[3,4]))==2 % Bottom left corner
                                Nd_c_1=NODE{N_I+N_L-1};
                                Nd_c_2=NODE{N};
                                Nd_c_3=NODE{N_I+N_L-1+N_H-1};
                                %
                                cell_pool_2_2_1=Nd_c_1{5};
                                Cell_pool_2_2_1=cell(1,length(cell_pool_2_2_1));
                                for s=1:length(cell_pool_2_2_1)
                                    C=Cell_pool_2_2_1{s};
                                    Cell_p=CELL{cell_pool_2_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);-(Y2-Y1)];
                                    Cell_pool_2_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2_2=Nd_c_2{5};
                                Cell_pool_2_2_2=cell(1,length(cell_pool_2_2_2));
                                for s=1:length(cell_pool_2_2_2)
                                    C=Cell_pool_2_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                    Cell_pool_2_2_2{s}=C;
                                end
                                %
                                cell_pool_2_2_3=Nd_c_3{5};
                                Cell_pool_2_2_3=cell(1,length(cell_pool_2_2_3));
                                for s=1:length(cell_pool_2_2_3)
                                    C=Cell_pool_2_2_3{s};
                                    Cell_p=CELL{cell_pool_2_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                    Cell_pool_2_2_3{s}=C;
                                end
                            elseif length(union(Onb_c,[4,1]))==2 % Top left corner
                                Nd_c_1=NODE{N_I+N_L-1+N_H-1};
                                Nd_c_2=NODE{N_I+N_L-1};
                                Nd_c_3=NODE{N_I+N_L-1+N_H-2+N_L};
                                %
                                cell_pool_2_2_1=Nd_c_1{5};
                                Cell_pool_2_2_1=cell(1,length(cell_pool_2_2_1));
                                for s=1:length(cell_pool_2_2_1)
                                    C=Cell_pool_2_2_1{s};
                                    Cell_p=CELL{cell_pool_2_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);(Y2-Y1)];
                                    Cell_pool_2_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2_2=Nd_c_2{5};
                                Cell_pool_2_2_2=cell(1,length(cell_pool_2_2_2));
                                for s=1:length(cell_pool_2_2_2)
                                    C=Cell_pool_2_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                    Cell_pool_2_2_2{s}=C;
                                end
                                %
                                cell_pool_2_2_3=Nd_c_3{5};
                                Cell_pool_2_2_3=cell(1,length(cell_pool_2_2_3));
                                for s=1:length(cell_pool_2_2_3)
                                    C=Cell_pool_2_2_3{s};
                                    Cell_p=CELL{cell_pool_2_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                    Cell_pool_2_2_3{s}=C;
                                end
                            else
                                error('logic error1');
                            end
                            %% Combine all pools
                            if length(union(cell_pool_2_2_1,union(cell_pool_2_2_2,cell_pool_2_2_3)))~=(length(cell_pool_2_2_1)+length(cell_pool_2_2_2)+length(cell_pool_2_2_3))
                                error('The three sub pools should be mutually exclusive!');
                            end
                            cell_pool_2_2=[cell_pool_2_2_1,cell_pool_2_2_2,cell_pool_2_2_3];
                            Cell_pool_2_2={Cell_pool_2_2_1{1:end},Cell_pool_2_2_2{1:end},Cell_pool_2_2_3{1:end}}; % It is important to use {1:end}, otherwise, the entire value will be on entry of the new cell
                            if isempty(intersect(cell_pool_2_1,cell_pool_2_2))
                                error('Logic error!');
                            end
                            cell_pool_2=union(cell_pool_2_1,cell_pool_2_2);
                            if length(cell_pool_2)>=(length(cell_pool_2_1)+length(cell_pool_2_2))
                                error('Logic error!');
                            end
                            cell_pool_2_temp=[cell_pool_2_1,cell_pool_2_2];
                            Cell_pool_2_temp={Cell_pool_2_1{1:end},Cell_pool_2_2{1:end}};
                            Cell_pool_2=cell(1,length(cell_pool_2));
                            for s=1:length(cell_pool_2)
                                for j=1:length(cell_pool_2_temp)
                                    if cell_pool_2(s)==cell_pool_2_temp(j)
                                        Cell_pool_2{s}=Cell_pool_2_temp{j};
                                        break;
                                    end
                                end
                            end
                            if ~isempty(intersect(cell_pool_1,cell_pool_2))
                                error('Logic error!');
                            end
                            cell_pool=[cell_pool_1,cell_pool_2];
                            Cell_pool={Cell_pool_1{1:end},Cell_pool_2{1:end}};
                            for s=1:length(cell_pool)
                                C_p=Cell_pool{s};
                                if cell_pool(s)~=C_p{1}
                                    error('Not match!');
                                end
                            end
                        elseif length(Onb1)==1 && length(Onb2)==1
                            %% First centroid pool
                            cell_pool_1=setxor(Cell_fix{1},union(Nd1{5},Nd2{5})); % The centroid candidates on the local side of the periodic boundaries
                            Cell_pool_1=cell(1,length(cell_pool_1)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                            for s=1:length(cell_pool_1)
                                C=Cell_pool_1{s};
                                Cell_p=CELL{cell_pool_1(s)};
                                C{1,1}=Cell_p{1};
                                C{2,1}=Cell_p{5};
                                Cell_pool_1{s}=C;
                            end
                            %% Second centroid pool
                            cell_pool_2=union(setxor(Nd1{5},Nd1{9}),setxor(Nd2{5},Nd2{9})); % The centroid candidates on the periodic side of the periodic boundaries
                            Cell_pool_2=cell(1,length(cell_pool_2)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                            if Onb1==1 && Onb2==1
                                for s=1:length(cell_pool_2)
                                    C=Cell_pool_2{s};
                                    Cell_p=CELL{cell_pool_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                    Cell_pool_2{s}=C;
                                end
                            elseif Onb1==2 && Onb2==2
                                for s=1:length(cell_pool_2)
                                    C=Cell_pool_2{s};
                                    Cell_p=CELL{cell_pool_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                    Cell_pool_2{s}=C;
                                end
                            elseif Onb1==3 && Onb2==3
                                for s=1:length(cell_pool_2)
                                    C=Cell_pool_2{s};
                                    Cell_p=CELL{cell_pool_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                    Cell_pool_2{s}=C;
                                end
                            elseif Onb1==4 && Onb2==4
                                for s=1:length(cell_pool_2)
                                    C=Cell_pool_2{s};
                                    Cell_p=CELL{cell_pool_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                    Cell_pool_2{s}=C;
                                end
                            else
                                error('Logic error1');
                            end
                            %% Combine pools
                            if ~isempty(intersect(cell_pool_1,cell_pool_2))
                                error('Logic error!');
                            end
                            cell_pool=[cell_pool_1,cell_pool_2];
                            Cell_pool={Cell_pool_1{1:end},Cell_pool_2{1:end}};
                            for s=1:length(cell_pool)
                                C_p=Cell_pool{s};
                                if cell_pool(s)~=C_p{1}
                                    error('Not match!');
                                end
                            end
                        else
                            error('Logic error!');
                        end
                    elseif (Face_target{10}~=0 && Face_target{11}==0) || (Face_target{10}==0 && Face_target{11}~=0) % one of two end nodes of target face are located on outer boundary
                        if Face_target{10}~=0
                            Onb=on_which_boundary(Nd1{3},X1,X2,Y1,Y2);
                            Nd_b=Nd1;
                        elseif Face_target{11}~=0
                            Onb=on_which_boundary(Nd2{3},X1,X2,Y1,Y2);
                            Nd_b=Nd2;
                        else
                            error('Logic error!');
                        end
                        if length(Onb)==2
                            %% First centroid pool
                            cell_pool_1=setxor(Cell_fix{1},union(Nd1{5},Nd2{5})); % The centroid candidates on the local side of the periodic boundaries
                            Cell_pool_1=cell(1,length(cell_pool_1)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                            for s=1:length(cell_pool_1)
                                C=Cell_pool_1{s};
                                Cell_p=CELL{cell_pool_1(s)};
                                C{1,1}=Cell_p{1};
                                C{2,1}=Cell_p{5};
                                Cell_pool_1{s}=C;
                            end
                            %% Second centroid pool
                            if length(union(Onb,[1,2]))==2 % Top right corner
                                Nd_c_1=NODE{N_I+N_L-1+N_H-2+N_L};
                                Nd_c_2=NODE{N_I+N_L-1+N_H-1};
                                Nd_c_3=NODE{N};
                                %
                                cell_pool_2_1=Nd_c_1{5};
                                Cell_pool_2_1=cell(1,length(cell_pool_2_1));
                                for s=1:length(cell_pool_2_1)
                                    C=Cell_pool_2_1{s};
                                    Cell_p=CELL{cell_pool_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);(Y2-Y1)];
                                    Cell_pool_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2=Nd_c_2{5};
                                Cell_pool_2_2=cell(1,length(cell_pool_2_2));
                                for s=1:length(cell_pool_2_2)
                                    C=Cell_pool_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                    Cell_pool_2_2{s}=C;
                                end
                                %
                                cell_pool_2_3=Nd_c_3{5};
                                Cell_pool_2_3=cell(1,length(cell_pool_2_3));
                                for s=1:length(cell_pool_2_3)
                                    C=Cell_pool_2_3{s};
                                    Cell_p=CELL{cell_pool_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                    Cell_pool_2_3{s}=C;
                                end
                            elseif length(union(Onb,[2,3]))==2 % Botom right corner
                                Nd_c_1=NODE{N};
                                Nd_c_2=NODE{N_I+N_L-1+N_H-2+N_L};
                                Nd_c_3=NODE{N_I+N_L-1};
                                %
                                cell_pool_2_1=Nd_c_1{5};
                                Cell_pool_2_1=cell(1,length(cell_pool_2_1));
                                for s=1:length(cell_pool_2_1)
                                    C=Cell_pool_2_1{s};
                                    Cell_p=CELL{cell_pool_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);-(Y2-Y1)];
                                    Cell_pool_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2=Nd_c_2{5};
                                Cell_pool_2_2=cell(1,length(cell_pool_2_2));
                                for s=1:length(cell_pool_2_2)
                                    C=Cell_pool_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[(X2-X1);0];
                                    Cell_pool_2_2{s}=C;
                                end
                                %
                                cell_pool_2_3=Nd_c_3{5};
                                Cell_pool_2_3=cell(1,length(cell_pool_2_3));
                                for s=1:length(cell_pool_2_3)
                                    C=Cell_pool_2_3{s};
                                    Cell_p=CELL{cell_pool_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                    Cell_pool_2_3{s}=C;
                                end
                            elseif length(union(Onb,[3,4]))==2 % Bottom left corner
                                Nd_c_1=NODE{N_I+N_L-1};
                                Nd_c_2=NODE{N};
                                Nd_c_3=NODE{N_I+N_L-1+N_H-1};
                                %
                                cell_pool_2_1=Nd_c_1{5};
                                Cell_pool_2_1=cell(1,length(cell_pool_2_1));
                                for s=1:length(cell_pool_2_1)
                                    C=Cell_pool_2_1{s};
                                    Cell_p=CELL{cell_pool_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);-(Y2-Y1)];
                                    Cell_pool_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2=Nd_c_2{5};
                                Cell_pool_2_2=cell(1,length(cell_pool_2_2));
                                for s=1:length(cell_pool_2_2)
                                    C=Cell_pool_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                    Cell_pool_2_2{s}=C;
                                end
                                %
                                cell_pool_2_3=Nd_c_3{5};
                                Cell_pool_2_3=cell(1,length(cell_pool_2_3));
                                for s=1:length(cell_pool_2_3)
                                    C=Cell_pool_2_3{s};
                                    Cell_p=CELL{cell_pool_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                    Cell_pool_2_3{s}=C;
                                end
                            elseif length(union(Onb,[4,1]))==2 % Top left corner
                                Nd_c_1=NODE{N_I+N_L-1+N_H-1};
                                Nd_c_2=NODE{N_I+N_L-1};
                                Nd_c_3=NODE{N_I+N_L-1+N_H-2+N_L};
                                %
                                cell_pool_2_1=Nd_c_1{5};
                                Cell_pool_2_1=cell(1,length(cell_pool_2_1));
                                for s=1:length(cell_pool_2_1)
                                    C=Cell_pool_2_1{s};
                                    Cell_p=CELL{cell_pool_2_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);(Y2-Y1)];
                                    Cell_pool_2_1{s}=C;
                                end
                                %
                                cell_pool_2_2=Nd_c_2{5};
                                Cell_pool_2_2=cell(1,length(cell_pool_2_2));
                                for s=1:length(cell_pool_2_2)
                                    C=Cell_pool_2_2{s};
                                    Cell_p=CELL{cell_pool_2_2(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                    Cell_pool_2_2{s}=C;
                                end
                                %
                                cell_pool_2_3=Nd_c_3{5};
                                Cell_pool_2_3=cell(1,length(cell_pool_2_3));
                                for s=1:length(cell_pool_2_3)
                                    C=Cell_pool_2_3{s};
                                    Cell_p=CELL{cell_pool_2_3(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                    Cell_pool_2_3{s}=C;
                                end
                            else
                                error('logic error1');
                            end
                            %% Combine all pools
                            if length(union(cell_pool_2_1,union(cell_pool_2_2,cell_pool_2_3)))~=(length(cell_pool_2_1)+length(cell_pool_2_2)+length(cell_pool_2_3))
                                error('The three sub pools should be mutually exclusive!');
                            end
                            cell_pool_2=[cell_pool_2_1,cell_pool_2_2,cell_pool_2_3];
                            Cell_pool_2={Cell_pool_2_1{1:end},Cell_pool_2_2{1:end},Cell_pool_2_3{1:end}}; % It is important to use {1:end}, otherwise, the entire value will be on entry of the new cell
                            if ~isempty(intersect(cell_pool_1,cell_pool_2))
                                error('Logic error!');
                            end
                            cell_pool=[cell_pool_1,cell_pool_2];
                            Cell_pool={Cell_pool_1{1:end},Cell_pool_2{1:end}};
                            for s=1:length(cell_pool)
                                C_p=Cell_pool{s};
                                if cell_pool(s)~=C_p{1}
                                    error('Not match!');
                                end
                            end
                        else
                            if Onb==0
                                error('This node should be on boundary!');
                            else
                                %% First centroid pool
                                cell_pool_1=setxor(Cell_fix{1},union(Nd1{5},Nd2{5})); % The centroid candidates on the local side of the periodic boundaries
                                Cell_pool_1=cell(1,length(cell_pool_1)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                                for s=1:length(cell_pool_1)
                                    C=Cell_pool_1{s};
                                    Cell_p=CELL{cell_pool_1(s)};
                                    C{1,1}=Cell_p{1};
                                    C{2,1}=Cell_p{5};
                                    Cell_pool_1{s}=C;
                                end
                                %% Second centroid pool
                                cell_pool_2=setxor(Nd_b{5},Nd_b{9}); % The centroid candidates on the periodic side of the periodic boundaries
                                Cell_pool_2=cell(1,length(cell_pool_2)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                                if Onb==1
                                    for s=1:length(cell_pool_2)
                                        C=Cell_pool_2{s};
                                        Cell_p=CELL{cell_pool_2(s)};
                                        C{1,1}=Cell_p{1};
                                        C{2,1}=Cell_p{5}+[0;(Y2-Y1)];
                                        Cell_pool_2{s}=C;
                                    end
                                elseif Onb==2
                                    for s=1:length(cell_pool_2)
                                        C=Cell_pool_2{s};
                                        Cell_p=CELL{cell_pool_2(s)};
                                        C{1,1}=Cell_p{1};
                                        C{2,1}=Cell_p{5}+[(X2-X1);0];
                                        Cell_pool_2{s}=C;
                                    end
                                elseif Onb==3
                                    for s=1:length(cell_pool_2)
                                        C=Cell_pool_2{s};
                                        Cell_p=CELL{cell_pool_2(s)};
                                        C{1,1}=Cell_p{1};
                                        C{2,1}=Cell_p{5}+[0;-(Y2-Y1)];
                                        Cell_pool_2{s}=C;
                                    end
                                elseif Onb==4
                                    for s=1:length(cell_pool_2)
                                        C=Cell_pool_2{s};
                                        Cell_p=CELL{cell_pool_2(s)};
                                        C{1,1}=Cell_p{1};
                                        C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                        Cell_pool_2{s}=C;
                                    end
                                else
                                    error('Logic error!');
                                end
                                %% Combine pools
                                if ~isempty(intersect(cell_pool_1,cell_pool_2))
                                    error('Logic error!');
                                end
                                cell_pool=[cell_pool_1,cell_pool_2];
                                Cell_pool={Cell_pool_1{1:end},Cell_pool_2{1:end}};
                                for s=1:length(cell_pool)
                                    C_p=Cell_pool{s};
                                    if cell_pool(s)~=C_p{1}
                                        error('Not match!');
                                    end
                                end
                            end
                        end
                    elseif Face_target{10}==0 && Face_target{11}==0 % Both end nodes of target face are interior, Copy FC{20}
                        cell_pool=setxor(Cell_fix{1},union(Nd1{5},Nd2{5})); % The centroid candidates on the local side of the periodic boundaries
                        Cell_pool=cell(1,length(cell_pool)); % Temperary cell structure, which contains only the cell order number and centroid coordinates
                        for s=1:length(cell_pool)
                            C=Cell_pool{s};
                            Cell_p=CELL{cell_pool(s)};
                            C{1,1}=Cell_p{1};
                            C{2,1}=Cell_p{5};
                            Cell_pool{s}=C;
                        end
                        for s=1:length(cell_pool)
                            C_p=Cell_pool{s};
                            if cell_pool(s)~=C_p{1}
                                error('Not match!');
                            end
                        end
                    else
                        error('Logic error!');
                    end
                    %% Find all possible paired centroids that could be used to form triangles
                    L=length(cell_pool);
                    cell_pool_count=1:1:L;
                    Pair_cell_union=zeros(2,L*(L-1)/2);
                    a=0;
                    for k=1:length(cell_pool_count)
                        if length(cell_pool_count)==2
                            a=a+1;
                            Pair_cell_union(:,a)=[cell_pool_count(1);cell_pool_count(2)];
                            break;
                        else
                            for n=2:length(cell_pool_count)
                                a=a+1;
                                Pair_cell_union(:,a)=[cell_pool_count(1);cell_pool_count(n)];
                            end
                            cell_pool_count=setxor(cell_pool_count(1),cell_pool_count);
                        end
                    end
                    % check
                    if a~=L*(L-1)/2
                        error('Logic error!');
                    end
                    %% Narrow down the previous possibility ruling out the triangles that don't circle the stencil point in the upstream cell
                    b=0;
                    in_tri_pair=0;
                    for k=1:a
                        Cell_pair_1=Cell_pool{Pair_cell_union(1,k)};
                        Cell_pair_2=Cell_pool{Pair_cell_union(2,k)};
                        C_1=Cell_pair_1{2};
                        C_2=Cell_pair_2{2};
                        if in_triangle(fc26(:,i),C_fix,C_1,C_2)
                            b=b+1;
                            in_tri_pair(b)=k;
                        end
                    end
                    if b==0
                        error('No triangle is found that contains the stencil point!');
                    end
                    %% Find the cell that has the shortest distance to the stencil point
                    Dis=zeros(1,b);
                    for k=1:b
                        Cell_pair_1=Cell_pool{Pair_cell_union(1,in_tri_pair(k))};
                        Cell_pair_2=Cell_pool{Pair_cell_union(2,in_tri_pair(k))};
                        Dis(1,k)=(dis(fc26(:,i),Cell_pair_1{2})+dis(fc26(:,i),Cell_pair_2{2}))/2;
                    end
                    Dis_min=min(Dis);
                    for k=1:b
                        if single(e+Dis_min)==single(e+Dis(1,k))
                            break;
                        end
                    end
                    if k==b
                        if Dis(1,k)~=Dis_min
                            error('Logic error!');
                        end
                    end
                    
                    Cell_pair_1=Cell_pool{Pair_cell_union(1,in_tri_pair(k))};
                    Cell_pair_2=Cell_pool{Pair_cell_union(2,in_tri_pair(k))};
                    fc29(:,i)=[0;Cell_fix{1};Cell_pair_1{1};Cell_pair_2{1}];
                    fc30(:,i)=[C_fix;Cell_pair_1{2};Cell_pair_2{2}];
                    %% Check
                    if ~in_triangle(fc26(:,i),fc30(1:2,i),fc30(3:4,i),fc30(5:6,i))
                        error('Logic error!');
                    end
                end
            end
            FC{29}=fc29;
            FC{30}=fc30;
            %% Final data filling
            FACE{l}=FC;
        end
    end
    % Check FC{29},FC{30}
    C_TL=[X1;Y2];
    C_TR=[X2;Y2];
    C_BR=[X2;Y1];
    C_BL=[X1;Y1];
    for l=1:O
        FC=FACE{l};
        fc24=FC{24};
        if fc24==2 || (fc24==3 || fc24==4)  % Face type 2, 3 and 4
            fc25=FC{25};
            fc26=FC{26};
            fc28=FC{28};
            fc29=FC{29};
            fc30=FC{30};
            for i=1:length(fc25)
                if fc29(1,i)~=0
                    error('Should be zero!');
                end
                
                if fc28(1,i)==4
                    if fc29(2,i)~=fc29(3,i) || (fc29(3,i)~=fc29(4,i) || fc29(4,i)~=fc29(2,i))
                        error('FC{29} is incorrect!');
                    end
                    if single(e+sum(fc30(1:2,i)-fc30(3:4,i)))~=single(e) || (single(e+sum(fc30(3:4,i)-fc30(5:6,i)))~=single(e) || single(e+sum(fc30(5:6,i)-fc30(1:2,i)))~=single(e))
                        error('FC{30} is incorrect!');
                    end
                else
                    if ~in_triangle(fc26(:,i),fc30(1:2,i),fc30(3:4,i),fc30(5:6,i))
                        error('The stencil point is not enclosed by found three centroids!');
                    end
                    Cell_fix=CELL{fc29(2,i)};
                    C_fix=Cell_fix{5};
                    C_fix_t=fc30(1:2,i);
                    if single(e+dis(C_fix,C_fix_t))>single(e)
                        error('The coordinate in FC{30} is incorrect!');
                    end
                    for j=1:2
                        Cell_cent=CELL{fc29(2+j,i)};
                        C_cent=Cell_cent{5};
                        C_cent_t=fc30(3+(j-1)*2:4+(j-1)*2,i);
                        if single(e+dis(C_cent,C_cent_t))>single(e)
%                             v=C_cent_t-fc26(:,i);
%                             [~, C_i_top]=intercept(fc26(:,i),v,C_TL,C_TR);
%                             [~, C_i_right]=intercept(fc26(:,i),v,C_TR,C_BR);
%                             [~, C_i_bottom]=intercept(fc26(:,i),v,C_BR,C_BL);
%                             [~, C_i_left]=intercept(fc26(:,i),v,C_BL,C_TL);
%                             on_top=on_edge(C_i_top,fc26(:,i),C_cent_t);
%                             on_right=on_edge(C_i_right,fc26(:,i),C_cent_t);
%                             on_bottom=on_edge(C_i_bottom,fc26(:,i),C_cent_t);
%                             on_left=on_edge(C_i_left,fc26(:,i),C_cent_t);
                            
                            v_top=C_TL-C_TR;
                            Ori_top=C_TR;
                            v_right=C_TR-C_BR;
                            Ori_right=C_BR;
                            v_bottom=C_BR-C_BL;
                            Ori_bottom=C_BL;
                            v_left=C_BL-C_TL;
                            Ori_left=C_TL;
                            [on_top, ~]=intercept(Ori_top,v_top,fc26(:,i),C_cent_t);
                            [on_right, ~]=intercept(Ori_right,v_right,fc26(:,i),C_cent_t);
                            [on_bottom, ~]=intercept(Ori_bottom,v_bottom,fc26(:,i),C_cent_t);
                            [on_left, ~]=intercept(Ori_left,v_left,fc26(:,i),C_cent_t);
%                             on_top=on_edge(C_i_top,C_TL,C_TR);
%                             on_right=on_edge(C_i_right,C_TR,C_BR);
%                             on_bottom=on_edge(C_i_bottom,C_BR,C_BL);
%                             on_left=on_edge(C_i_left,C_BL,C_TL);
                            if (on_top+on_right+on_bottom+on_left)~=1
                                error('There is only one boundary being across!');
                            end
                            if on_top
                                if single(e+sum(C_cent+[0;(Y2-Y1)]-C_cent_t))~=single(e)
                                    error('FC{29} and/or FC{30} contains incorrect info!');
                                end
                            elseif on_right
                                if single(e+sum(C_cent+[(X2-X1);0]-C_cent_t))~=single(e)
                                    error('FC{29} and/or FC{30} contains incorrect info!');
                                end
                            elseif on_bottom
                                if single(e+sum(C_cent+[0;-(Y2-Y1)]-C_cent_t))~=single(e)
                                    error('FC{29} and/or FC{30} contains incorrect info!');
                                end
                            elseif on_left
                                if single(e+sum(C_cent+[-(X2-X1);0]-C_cent_t))~=single(e)
                                    error('FC{29} and/or FC{30} contains incorrect info!');
                                end
                            else
                                error('Logic error!');
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% Final check for IRT mesh
    if FM==0
        %% FC{16} FC{17}
        for l=(N_L-1)*(N_H-1)*4+1:N_I
            FC=FACE{l};
            fc16=FC{16};
            fc17=FC{17};
            for i=1:length(fc16)
                if fc16(1,i)~=0
                    P=CELL{fc16(1,i)};
                    if single(e1+sum(P{5}-fc17(:,i)))~=single(e1)
                        error('Incorrect FC{16} and/or FC{17} for IRT mesh!');
                    end
                end
            end
        end
        %% FC{20}
%         L=0;
%         L_counter=0;
%         for l=1:O
%             FC=FACE{l};
%             fc16=FC{16};
%             fc20=FC{20};
%             for i=1:length(fc16)
%                 if fc16(1,i)~=0
%                     if fc20(1,i)==1
%                         ND=NODE{fc20(2,i)};
%                         P1=CELL{fc20(3,i)};
%                         P2=CELL{fc20(4,i)};
%                         L_counter=L_counter+1;
%                         L(L_counter)=double(dis(ND{3},P1{5})+dis(ND{3},P2{5})+dis(P1{5},P2{5})+e1);
%                     end
%                 end
%             end
%         end
%         if length(unique(L))~=2
%             length(unique(L))
%             unique(L)
%             error('The stencil point is the on the edge of the enclosing triangle!');
%         end
        %% FC{25} FC{26}
        for l=1:O
            FC=FACE{l};
            fc25=FC{25};
            fc26=FC{26};
            if FC{24}==2
                for i=1:length(fc25)
                    P=CELL{fc25(1,i)};
                    if single(e1+sum(P{5}-fc26(:,i)))~=single(e1)
                        error('Incorrect FC{25} and/or FC{26} for IRT mesh!');
                    end
                end
            end
        end
    end
end

%% FC{31}~FC{37}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3 || FC{24}==4
        break;
    end
end
if ((isempty(FC{31}) || isempty(FC{32})) || (isempty(FC{33}) || isempty(FC{34}))) ...
        || ((isempty(FC{35}) || isempty(FC{36})) || isempty(FC{37}))
    for l=1:O
        FC=FACE{l};
        if FC{24}==2 || (FC{24}==3 || FC{24}==4)
            fc31=FC{16};
            fc32=FC{17};
            fc33=FC{22};
            fc34=FC{18};
            fc35=FC{19};
            fc36=FC{20};
            fc37=zeros(6,length(fc31));
            for i=1:length(fc31)
                if fc31(1,i)~=0
                    if fc36(1,i)==0
                        Cell1=CELL{fc36(2,i)};
                        Cell2=CELL{fc36(3,i)};
                        Cell3=CELL{fc36(4,i)};
                        P1=Cell1{5};
                        P2=Cell2{5};
                        P3=Cell3{5};
                    elseif fc36(1,i)==1
                        Node1=NODE{fc36(2,i)};
                        Cell2=CELL{fc36(3,i)};
                        Cell3=CELL{fc36(4,i)};
                        P1=Node1{3};
                        P2=Cell2{5};
                        P3=Cell3{5};
                    elseif fc36(1,i)==2
                        Node1=NODE{fc36(2,i)};
                        Node2=NODE{fc36(3,i)};
                        Cell3=CELL{fc36(4,i)};
                        P1=Node1{3};
                        P2=Node2{3};
                        P3=Cell3{5};
                    else
                        error('The stencil point cannot be enclosed by all boundary nodes!');
                    end
                    fc37(:,i)=[P1;P2;P3];
                end
            end
            fc23=FC{23};
            fc25=FC{25};
            fc26=FC{26};
            fc27=FC{27};
            fc28=FC{28};
            fc29=FC{29};
            fc30=FC{30};
            for i=1:length(fc23)
                if fc23(1,i)==2 || fc23(1,i)==4 % check whether the stencil point across left and right boundary, if yes, replace everything for that stancil point
                    fc31(1,i)=fc25(1,i);
                    fc32(:,i)=fc26(:,i);
                    fc33(1,i)=fc27(1,i);
                    fc34(:,i)=zeros(length(fc34(:,i)),1);
                    fc35(1,i)=fc28(1,i);
                    % Determine fc36 and fc37
                    Cell=CELL{fc31(1,i)};
                    if fc35(1,i)==4
                        fc36(:,i)=fc29(:,i);
                        fc37(:,i)=fc30(:,i);
                    else
                        if fc35(1,i)==1
                            ND1=NODE{Cell{7}};
                            ND2=NODE{Cell{8}};
                        elseif fc35(1,i)==2
                            ND1=NODE{Cell{8}};
                            ND2=NODE{Cell{9}};
                        elseif fc35(1,i)==3
                            ND1=NODE{Cell{9}};
                            ND2=NODE{Cell{7}};
                        else
                            error('Logic error!');
                        end
                        if (((ND1{1}==N_I+N_L-1 || ND1{1}==N_I+N_L-1+N_H-1) || (ND1{1}==N_I+N_L-1+N_H-1+N_L-1 || ND1{1}==N)) && ND2{1}<=N_I) || ...
                                (((ND2{1}==N_I+N_L-1 || ND2{1}==N_I+N_L-1+N_H-1) || (ND2{1}==N_I+N_L-1+N_H-1+N_L-1 || ND2{1}==N)) && ND1{1}<=N_I)
                            if ((ND1{1}==N_I+N_L-1 || ND1{1}==N_I+N_L-1+N_H-1) || (ND1{1}==N_I+N_L-1+N_H-1+N_L-1 || ND1{1}==N)) && ND2{1}<=N_I
                                ND=ND1;
                                C_nd=ND1{3};
                            else
                                ND=ND2;
                                C_nd=ND2{3};
                            end
                            C_fixed=Cell{5};
                            Cell_third_pool=union(ND1{5},ND2{5});
                            Cell_third_pool=setxor(Cell{1},Cell_third_pool);
                            L=length(Cell_third_pool);
                            a=0;
                            in_tri=0;
                            for k=1:L
                                Cell_third=CELL{Cell_third_pool(k)};
                                C_third_cell=Cell_third{5};
                                if in_triangle(fc32(:,i),C_fixed,C_nd,C_third_cell)
                                    a=a+1;
                                    in_tri(a)=k;
                                end
                            end
                            if a==0
                                error('Enclosing triangle is not found!');
                            end
                            % Find the cell that has the shortest distance to the stencil point
                            Dis=zeros(1,a);
                            for k=1:a
                                Cell_third=CELL{Cell_third_pool(in_tri(k))};
                                C_third_cell=Cell_third{5};
                                Dis(1,k)=dis(fc32(:,i),C_third_cell);
                            end
                            Dis_min=min(Dis);
                            for k=1:a
                                if double(e+Dis_min)==double(e+Dis(1,k))
                                    break;
                                end
                            end
                            if k==a
                                if Dis(1,k)~=Dis_min
                                    error('Logic error!');
                                end
                            end
                            Third_cell_found=Cell_third_pool(in_tri(k));
                            % Check
                            Cell_third=CELL{Third_cell_found};
                            C_third_cell=Cell_third{5};
                            if ~in_triangle(fc32(:,i),C_fixed,C_nd,C_third_cell)
                                error('Enclosing triangle is not found!');
                            end
                            fc36(:,i)=[1;ND{1};fc31(1,i);Third_cell_found];
                            fc37(:,i)=[C_nd;C_fixed;C_third_cell];
                        else % The stencil point is not close to the four corners
                            fc36(:,i)=fc29(:,i);
                            fc37(:,i)=fc30(:,i);
                        end
                    end
                else
                    if fc36(1,i)==1
                        ND=NODE{fc36(2,i)};
                        Onb=on_which_boundary(ND{3},X1,X2,Y1,Y2);
                        if length(Onb)==1
                            if Onb==2 || Onb==4
                                fc36(:,i)=fc29(:,i);
                                fc37(:,i)=fc30(:,i);
                            end
                        end
                    elseif fc36(1,i)==2
                        for s=1:2
                            ND=NODE{fc36(1+s,i)};
                            Onb=on_which_boundary(ND{3},X1,X2,Y1,Y2);
                            if length(Onb)==1
                                if Onb==2 || Onb==4
                                    fc36(:,i)=fc29(:,i);
                                    fc37(:,i)=fc30(:,i);
                                    break;
                                end
                            end
                        end
                    else
                        if fc36(1,i)~=0
                            error('Logic error!')
                        end
                    end
                end
                % Check FC{37} to make sure the three enclosing points
                % don't across top and bottom boundary
                y1=fc37(2,i);
                y2=fc37(4,i);
                y3=fc37(6,i);
                if abs(y1-y2)/(Y2-Y1)>0.5 || abs(y2-y3)/(Y2-Y1)>0.5 || abs(y3-y1)/(Y2-Y1)>0.5
                    error('Some of the enclosing points cross the top or bottom boundary!');
                end
            end
            % fill
            FC{31}=fc31;
            FC{32}=fc32;
            FC{33}=fc33;
            FC{34}=fc34;
            FC{35}=fc35;
            FC{36}=fc36;
            FC{37}=fc37;
            
            FACE{l}=FC;
        end
    end
    %% Check
    for l=1:O
        FC=FACE{l};
        if FC{24}==2 || (FC{24}==3 || FC{24}==4)
            fc31=FC{31};
            fc32=FC{32};
            fc33=FC{33};
            fc34=FC{34};
            fc35=FC{35};
            fc36=FC{36};
            fc37=FC{37};
            %% check fc31 and fc32
            for i=1:length(fc31)
                if fc31(1,i)~=0
                    C=CELL{fc31(1,i)};
                    C_1=C{13};
                    C_2=C{14};
                    C_3=C{15};
                    if ~in_triangle(fc32(:,i),C_1,C_2,C_3)
                        error('FC{31} and/or FC{32} contains false info');
                    end
                end
            end
            %% check fc32 and fc33
            if single(e1+dis(fc32(:,2),fc32(:,3)))~=single(e1+fc33(1,2)+fc33(1,3))
                if single(e1+dis(fc32(:,2)-[(X2-X1);0],fc32(:,3)))~=single(e1+fc33(1,2)+fc33(1,3))
                    if single(e1+dis(fc32(:,2)+[(X2-X1);0],fc32(:,3)))~=single(e1+fc33(1,2)+fc33(1,3))
                        error('FC{32} and/or FC{33} contains false info');
                    end
                end
            end
            if single(e1+dis(fc32(:,1),fc32(:,4)))~=single(e1+fc33(1,1)+fc33(1,4))
                if single(e1+dis(fc32(:,1)-[(X2-X1);0],fc32(:,4)))~=single(e1+fc33(1,1)+fc33(1,4))
                    if single(e1+dis(fc32(:,1)+[(X2-X1);0],fc32(:,4)))~=single(e1+fc33(1,1)+fc33(1,4))
                        error('FC{32} and/or FC{33} contains false info');
                    end
                end
            end
            %% fc34
            for i=1:length(fc31)
                if sum(fc34(:,i))~=0
                    Nd1=NODE{fc34(1,i)};
                    Nd2=NODE{fc34(2,i)};
                    Intercept=(1-fc34(3,i))*Nd1{3}+fc34(3,i)*Nd2{3};
                    if i==1
                        S1=fc32(:,1);
                        S2=fc32(:,2);
                    elseif i==2
                        S1=fc32(:,2);
                        S2=fc32(:,3);
                    elseif i==3
                        S1=fc32(:,3);
                        S2=fc32(:,2);
                    elseif i==4
                        S1=fc32(:,3);
                        S2=fc32(:,4);
                    else
                        error('Logic error!');
                    end
                    if ~on_edge(Intercept,S1,S2)
                        if ~on_edge(single(Intercept),single(S1),single(S2))
                            error('FC{34} contains false info');
                        end
                    end
                end
            end
            %% fc35, fc36 and fc37
            for i=1:length(fc31)
                if fc31(1,i)~=0
                    if fc35(1,i)==4
                        if fc36(1,i)~=0 || (fc36(2,i)~=fc36(3,i) || fc36(3,i)~=fc36(4,i))
                            error('FC{35} and/or FC{36} contains false info');
                        end
                        if single(e1+sum(fc37(1:2,i)-fc37(3:4,i)))~=single(e1) || single(e1+sum(fc37(3:4,i)-fc37(5:6,i)))~=single(e1)
                            error('FC{37} contains false info');
                        end
                    else
                        if ~in_triangle(fc32(:,i),fc37(1:2,i),fc37(3:4,i),fc37(5:6,i))
                            error('FC{37} contains false info');
                        end
                    end
                end
            end
            %% fc35, fc36 and fc37
            for i=1:length(fc31)
                if fc31(1,i)~=0
                    if fc35(1,i)~=4
                        if fc36(1,i)==0
                            Cell1=CELL{fc36(2,i)};
                            Cell2=CELL{fc36(3,i)};
                            Cell3=CELL{fc36(4,i)};
                            P1=Cell1{5};
                            P2=Cell2{5};
                            P3=Cell3{5};
                            Dis1=dis(P1,P2);
                            Dis2=dis(P2,P3);
                            Dis3=dis(P3,P1);
                            if Dis1>(X2-X1)/2 && Dis2>(X2-X1)/2 && Dis3<(X2-X1)/2
                                if P2(1,1)<(P1(1,1)+P3(1,1))/2
                                    P2(1,1)=P2(1,1)+(X2-X1);
                                else
                                    P2(1,1)=P2(1,1)-(X2-X1);
                                end
                            elseif Dis2>(X2-X1)/2 && Dis3>(X2-X1)/2 && Dis1<(X2-X1)/2
                                if P3(1,1)<(P1(1,1)+P2(1,1))/2
                                    P3(1,1)=P3(1,1)+(X2-X1);
                                else
                                    P3(1,1)=P3(1,1)-(X2-X1);
                                end
                            elseif Dis3>(X2-X1)/2 && Dis1>(X2-X1)/2 && Dis2<(X2-X1)/2
                                if P1(1,1)<(P2(1,1)+P3(1,1))/2
                                    P1(1,1)=P1(1,1)+(X2-X1);
                                else
                                    P1(1,1)=P1(1,1)-(X2-X1);
                                end
                            elseif Dis1<(X2-X1)/2 && Dis2<(X2-X1)/2 && Dis3<(X2-X1)/2
                                ;
                            else
                                error('Logic error!');
                            end
                        elseif fc36(1,i)==1
                            Node1=NODE{fc36(2,i)};
                            Cell2=CELL{fc36(3,i)};
                            Cell3=CELL{fc36(4,i)};
                            P1=Node1{3};
                            P2=Cell2{5};
                            P3=Cell3{5};
                        elseif fc36(1,i)==2
                            Node1=NODE{fc36(2,i)};
                            Node2=NODE{fc36(3,i)};
                            Cell3=CELL{fc36(4,i)};
                            P1=Node1{3};
                            P2=Node2{3};
                            P3=Cell3{5};
                        else
                            error('The stencil point cannot be enclosed by all boundary nodes!');
                        end
                        if ~in_triangle(fc32(:,i),P1,P2,P3)
                            if ~in_triangle(fc32(:,i)+[(X2-X1);0],P1,P2,P3)
                                if ~in_triangle(fc32(:,i)-[(X2-X1);0],P1,P2,P3)
                                    error('FC{36} and/or FC{37} contains false info');
                                end
                            end
                        end
                    end
                end
            end
            
        end
    end
end

%% FC{38}~FC{44}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3 || FC{24}==4
        break;
    end
end
if ((isempty(FC{38}) || isempty(FC{39})) || (isempty(FC{40}) || isempty(FC{41}))) ...
        || ((isempty(FC{42}) || isempty(FC{43})) || isempty(FC{44}))
    for l=1:O
        FC=FACE{l};
        if FC{24}==2 || (FC{24}==3 || FC{24}==4)
            fc38=FC{16};
            fc39=FC{17};
            fc40=FC{22};
            fc41=FC{18};
            fc42=FC{19};
            fc43=FC{20};
            fc44=zeros(6,length(fc38));
            for i=1:length(fc38)
                if fc38(1,i)~=0
                    if fc43(1,i)==0
                        Cell1=CELL{fc43(2,i)};
                        Cell2=CELL{fc43(3,i)};
                        Cell3=CELL{fc43(4,i)};
                        P1=Cell1{5};
                        P2=Cell2{5};
                        P3=Cell3{5};
                    elseif fc43(1,i)==1
                        Node1=NODE{fc43(2,i)};
                        Cell2=CELL{fc43(3,i)};
                        Cell3=CELL{fc43(4,i)};
                        P1=Node1{3};
                        P2=Cell2{5};
                        P3=Cell3{5};
                    elseif fc43(1,i)==2
                        Node1=NODE{fc43(2,i)};
                        Node2=NODE{fc43(3,i)};
                        Cell3=CELL{fc43(4,i)};
                        P1=Node1{3};
                        P2=Node2{3};
                        P3=Cell3{5};
                    else
                        error('The stencil point cannot be enclosed by all boundary nodes!');
                    end
                    fc44(:,i)=[P1;P2;P3];
                end
            end
            fc23=FC{23};
            fc25=FC{25};
            fc26=FC{26};
            fc27=FC{27};
            fc28=FC{28};
            fc29=FC{29};
            fc30=FC{30};
            for i=1:length(fc23)
                if fc23(1,i)==1 || fc23(1,i)==3 % check whether the stencil point across top and bottom boundary, if yes, replace everything for that stancil point
                    fc38(1,i)=fc25(1,i);
                    fc39(:,i)=fc26(:,i);
                    fc40(1,i)=fc27(1,i);
                    fc41(:,i)=zeros(length(fc41(:,i)),1);
                    fc42(1,i)=fc28(1,i);
                    % Determine fc43 and fc44
                    Cell=CELL{fc38(1,i)};
                    if fc42(1,i)==4
                        fc43(:,i)=fc29(:,i);
                        fc44(:,i)=fc30(:,i);
                    else
                        if fc42(1,i)==1
                            ND1=NODE{Cell{7}};
                            ND2=NODE{Cell{8}};
                        elseif fc42(1,i)==2
                            ND1=NODE{Cell{8}};
                            ND2=NODE{Cell{9}};
                        elseif fc42(1,i)==3
                            ND1=NODE{Cell{9}};
                            ND2=NODE{Cell{7}};
                        else
                            error('Logic error!');
                        end
                        if (((ND1{1}==N_I+N_L-1 || ND1{1}==N_I+N_L-1+N_H-1) || (ND1{1}==N_I+N_L-1+N_H-1+N_L-1 || ND1{1}==N)) && ND2{1}<=N_I) || ...
                                (((ND2{1}==N_I+N_L-1 || ND2{1}==N_I+N_L-1+N_H-1) || (ND2{1}==N_I+N_L-1+N_H-1+N_L-1 || ND2{1}==N)) && ND1{1}<=N_I)
                            if ((ND1{1}==N_I+N_L-1 || ND1{1}==N_I+N_L-1+N_H-1) || (ND1{1}==N_I+N_L-1+N_H-1+N_L-1 || ND1{1}==N)) && ND2{1}<=N_I
                                ND=ND1;
                                C_nd=ND1{3};
                            else
                                ND=ND2;
                                C_nd=ND2{3};
                            end
                            C_fixed=Cell{5};
                            Cell_third_pool=union(ND1{5},ND2{5});
                            Cell_third_pool=setxor(Cell{1},Cell_third_pool);
                            L=length(Cell_third_pool);
                            a=0;
                            in_tri=0;
                            for k=1:L
                                Cell_third=CELL{Cell_third_pool(k)};
                                C_third_cell=Cell_third{5};
                                if in_triangle(fc39(:,i),C_fixed,C_nd,C_third_cell)
                                    a=a+1;
                                    in_tri(a)=k;
                                end
                            end
                            if a==0
                                error('Enclosing triangle is not found!');
                            end
                            % Find the cell that has the shortest distance to the stencil point
                            Dis=zeros(1,a);
                            for k=1:a
                                Cell_third=CELL{Cell_third_pool(in_tri(k))};
                                C_third_cell=Cell_third{5};
                                Dis(1,k)=dis(fc39(:,i),C_third_cell);
                            end
                            Dis_min=min(Dis);
                            for k=1:a
                                if double(e+Dis_min)==double(e+Dis(1,k))
                                    break;
                                end
                            end
                            if k==a
                                if Dis(1,k)~=Dis_min
                                    error('Logic error!');
                                end
                            end
                            Third_cell_found=Cell_third_pool(in_tri(k));
                            % Check
                            Cell_third=CELL{Third_cell_found};
                            C_third_cell=Cell_third{5};
                            if ~in_triangle(fc39(:,i),C_fixed,C_nd,C_third_cell)
                                error('Enclosing triangle is not found!');
                            end
                            fc43(:,i)=[1;ND{1};fc38(1,i);Third_cell_found];
                            fc44(:,i)=[C_nd;C_fixed;C_third_cell];
                        else % The stencil point is not close to the four corners
                            fc43(:,i)=fc29(:,i);
                            fc44(:,i)=fc30(:,i);
                        end
                    end
                else
                    if fc43(1,i)==1
                        ND=NODE{fc43(2,i)};
                        Onb=on_which_boundary(ND{3},X1,X2,Y1,Y2);
                        if length(Onb)==1
                            if Onb==1 || Onb==3
                                fc43(:,i)=fc29(:,i);
                                fc44(:,i)=fc30(:,i);
                            end
                        end
                    elseif fc43(1,i)==2
                        for s=1:2
                            ND=NODE{fc43(1+s,i)};
                            Onb=on_which_boundary(ND{3},X1,X2,Y1,Y2);
                            if length(Onb)==1
                                if Onb==1 || Onb==3
                                    fc43(:,i)=fc29(:,i);
                                    fc44(:,i)=fc30(:,i);
                                    break;
                                end
                            end
                        end
                    else
                        if fc43(1,i)~=0
                            error('Logic error!')
                        end
                    end
                end
                % Check FC{37} to make sure the three enclosing points
                % don't across top and bottom boundary
                x1=fc44(1,i);
                x2=fc44(3,i);
                x3=fc44(5,i);
                if abs(x1-x2)/(X2-X1)>0.5 || abs(x2-x3)/(X2-X1)>0.5 || abs(x3-x1)/(X2-X1)>0.5
                    error('Some of the enclosing points cross the left or right boundary!');
                end
            end
            % fill
            FC{38}=fc38;
            FC{39}=fc39;
            FC{40}=fc40;
            FC{41}=fc41;
            FC{42}=fc42;
            FC{43}=fc43;
            FC{44}=fc44;
            
            FACE{l}=FC;
        end
    end
    %% Check
    for l=1:O
        FC=FACE{l};
        if FC{24}==2 || (FC{24}==3 || FC{24}==4)
            fc38=FC{38};
            fc39=FC{39};
            fc40=FC{40};
            fc41=FC{41};
            fc42=FC{42};
            fc43=FC{43};
            fc44=FC{44};
            %% check fc38 and fc39
            for i=1:length(fc38)
                if fc38(1,i)~=0
                    C=CELL{fc38(1,i)};
                    C_1=C{13};
                    C_2=C{14};
                    C_3=C{15};
                    if ~in_triangle(fc39(:,i),C_1,C_2,C_3)
                        error('FC{38} and/or FC{39} contains false info');
                    end
                end
            end
            %% check fc39 and fc40
            if single(e1+dis(fc39(:,2),fc39(:,3)))~=single(e1+fc40(1,2)+fc40(1,3))
                if single(e1+dis(fc39(:,2)-[0;(Y2-Y1)],fc39(:,3)))~=single(e1+fc40(1,2)+fc40(1,3))
                    if single(e1+dis(fc39(:,2)+[0;(Y2-Y1)],fc39(:,3)))~=single(e1+fc40(1,2)+fc40(1,3))
                        error('FC{39} and/or FC{40} contains false info');
                    end
                end
            end
            if single(e1+dis(fc39(:,1),fc39(:,4)))~=single(e1+fc40(1,1)+fc40(1,4))
                if single(e1+dis(fc39(:,1)-[0;(Y2-Y1)],fc39(:,4)))~=single(e1+fc40(1,1)+fc40(1,4))
                    if single(e1+dis(fc39(:,1)+[0;(Y2-Y1)],fc39(:,4)))~=single(e1+fc40(1,1)+fc40(1,4))
                        error('FC{39} and/or FC{40} contains false info');
                    end
                end
            end
            %% fc41
            for i=1:length(fc38)
                if sum(fc41(:,i))~=0
                    Nd1=NODE{fc41(1,i)};
                    Nd2=NODE{fc41(2,i)};
                    Intercept=(1-fc41(3,i))*Nd1{3}+fc41(3,i)*Nd2{3};
                    if i==1
                        S1=fc39(:,1);
                        S2=fc39(:,2);
                    elseif i==2
                        S1=fc39(:,2);
                        S2=fc39(:,3);
                    elseif i==3
                        S1=fc39(:,3);
                        S2=fc39(:,2);
                    elseif i==4
                        S1=fc39(:,3);
                        S2=fc39(:,4);
                    else
                        error('Logic error!');
                    end
                    if ~on_edge(Intercept,S1,S2)
                        if ~on_edge(single(Intercept),single(S1),single(S2))
                            error('FC{34} contains false info');
                        end
                    end
                end
            end
            %% fc42, fc43 and fc44
            for i=1:length(fc38)
                if fc38(1,i)~=0
                    if fc42(1,i)==4
                        if fc43(1,i)~=0 || (fc43(2,i)~=fc43(3,i) || fc43(3,i)~=fc43(4,i))
                            error('FC{42} and/or FC43} contains false info');
                        end
                        if single(e1+sum(fc44(1:2,i)-fc44(3:4,i)))~=single(e1) || single(e1+sum(fc44(3:4,i)-fc44(5:6,i)))~=single(e1)
                            error('FC{44} contains false info');
                        end
                    else
                        if ~in_triangle(fc39(:,i),fc44(1:2,i),fc44(3:4,i),fc44(5:6,i))
                            error('FC{44} contains false info');
                        end
                    end
                end
            end
            %% fc35, fc36 and fc37
            for i=1:length(fc38)
                if fc38(1,i)~=0
                    if fc42(1,i)~=4
                        if fc43(1,i)==0
                            Cell1=CELL{fc43(2,i)};
                            Cell2=CELL{fc43(3,i)};
                            Cell3=CELL{fc43(4,i)};
                            P1=Cell1{5};
                            P2=Cell2{5};
                            P3=Cell3{5};
                            Dis1=dis(P1,P2);
                            Dis2=dis(P2,P3);
                            Dis3=dis(P3,P1);
                            if Dis1>(Y2-Y1)/2 && Dis2>(Y2-Y1)/2 && Dis3<(Y2-Y1)/2
                                if P2(2,1)<(P1(2,1)+P3(2,1))/2
                                    P2(2,1)=P2(2,1)+(Y2-Y1);
                                else
                                    P2(2,1)=P2(2,1)-(Y2-Y1);
                                end
                            elseif Dis2>(Y2-Y1)/2 && Dis3>(Y2-Y1)/2 && Dis1<(Y2-Y1)/2
                                if P3(2,1)<(P1(2,1)+P2(2,1))/2
                                    P3(2,1)=P3(2,1)+(Y2-Y1);
                                else
                                    P3(2,1)=P3(2,1)-(Y2-Y1);
                                end
                            elseif Dis3>(Y2-Y1)/2 && Dis1>(Y2-Y1)/2 && Dis2<(Y2-Y1)/2
                                if P1(2,1)<(P2(2,1)+P3(2,1))/2
                                    P1(2,1)=P1(2,1)+(Y2-Y1);
                                else
                                    P1(2,1)=P1(2,1)-(Y2-Y1);
                                end
                            elseif Dis1<(Y2-Y1)/2 && Dis2<(Y2-Y1)/2 && Dis3<(Y2-Y1)/2
                                ;
                            else
                                error('Logic error!');
                            end
                        elseif fc43(1,i)==1
                            Node1=NODE{fc43(2,i)};
                            Cell2=CELL{fc43(3,i)};
                            Cell3=CELL{fc43(4,i)};
                            P1=Node1{3};
                            P2=Cell2{5};
                            P3=Cell3{5};
                        elseif fc43(1,i)==2
                            Node1=NODE{fc43(2,i)};
                            Node2=NODE{fc43(3,i)};
                            Cell3=CELL{fc43(4,i)};
                            P1=Node1{3};
                            P2=Node2{3};
                            P3=Cell3{5};
                        else
                            error('The stencil point cannot be enclosed by all boundary nodes!');
                        end
                        if ~in_triangle(fc39(:,i),P1,P2,P3)
                            if ~in_triangle(fc39(:,i)+[0;(Y2-Y1)],P1,P2,P3)
                                if ~in_triangle(fc39(:,i)-[0;(Y2-Y1)],P1,P2,P3)
                                    error('FC{43} and/or FC{44} contains false info');
                                end
                            end
                        end
                    end
                end
            end
            
        end
    end
end


%% FC{45}, created, to be updated according to a seperate function with a macro flag during the run
FC=FACE{1};
if isempty(FC{45})
    for l=1:O
        FC=FACE{l};
        fc16=FC{16};
        fc45=ones(1,length(fc16));
        FC{45}=fc45;
        FACE{l}=FC;
    end
end

%%%%%%%%%%%%%%%%%%%%% From now, only lattice-dependent info %%%%%%%%%%%%%%%%
%% FC{46}, FC{47}
% FC{46}
for l=1:O
    FC=FACE{l};
    fc16=FC{16};
    fc22=FC{22};
    L_V_1=FC{4}*V1;
    
    %% initialization
    S=cell(1,30);
    s1=zeros(1,q1,'single');
    s2=zeros(1,q1,'single');
    s3=zeros(1,q1,'single');
    
    s4=zeros(1,q1);
    s5=zeros(1,q1);
    s6=zeros(1,q1);
    
    s7=zeros(1,q1);
    s8=zeros(1,q1);
    
    s9=zeros(1,q1);
    s10=zeros(1,q1);
    
    s11=zeros(1,q1);
    s12=zeros(1,q1);
    
    s13=zeros(1,q1);
    s14=zeros(1,q1);
    
    s15=zeros(1,q1)';
    
    s16=zeros(1,q1)';
    s17=zeros(1,q1)';
    
    s18=zeros(1,q1)';
    s19=zeros(1,q1)';
    s20=zeros(1,q1)';
    %% S{1}, S{2} and S{3}
    for k=1:q1
        if L_V_1(k)>=0
            s1(1,k)=fc16(1,2);
            s2(1,k)=fc16(1,3);
            s3(1,k)=fc16(1,4);
        else
            s1(1,k)=fc16(1,3);
            s2(1,k)=fc16(1,2);
            s3(1,k)=fc16(1,1);
        end
    end
    S{1}=s1;
    S{2}=s2;
    S{3}=s3;
    
    %% S{4}, S{5} and S{6}
    for k=1:q1
        if L_V_1(k)>=0
            s4(1,k)=fc22(1,2);
            s5(1,k)=fc22(1,3);
            s6(1,k)=fc22(1,4);
        else
            s4(1,k)=fc22(1,3);
            s5(1,k)=fc22(1,2);
            s6(1,k)=fc22(1,1);
        end
    end
    S{4}=s4;
    S{5}=s5;
    S{6}=s6;
    
    %% S{7~8}
    for k=1:q1
        s7(1,k)=fc22(1,2)+fc22(1,3);
        if L_V_1(k)>=0
            s8(1,k)=fc22(1,3)/s7(1,k);
        else
            s8(1,k)=fc22(1,2)/s7(1,k);
        end
    end
    S{7}=1./s7;
    S{8}=s8;
    
    %% S{9~10}
    for k=1:q1
        if L_V_1(k)>=0
            s9(1,k)=fc22(1,4)-fc22(1,3);
            s10(1,k)=fc22(1,3)/s9(1,k);
        else
            s9(1,k)=fc22(1,1)-fc22(1,2);
            s10(1,k)=fc22(1,2)/s9(1,k);
        end
    end
    S{9}=1./s9;
    S{10}=s10;
    
    %% S{11~12}
    for k=1:q1
        if L_V_1(k)>=0
            s11(1,k)=fc22(1,2)+fc22(1,4);
            s12(1,k)=fc22(1,3)/s11(1,k);
        else
            s11(1,k)=fc22(1,1)+fc22(1,3);
            s12(1,k)=fc22(1,2)/s11(1,k);
        end
    end
    S{11}=1./s11;
    S{12}=s12;
    
    %% S{13~14}
    for k=1:q1
        if L_V_1(k)>=0
            s13(1,k)=fc22(1,3);
        else
            s13(1,k)=fc22(1,2);
        end
    end
    S{13}=1./s13;
    S{14}=s14+1;
    
    %% S{15~20}
    for k=1:q1
        if L_V_1(k)>=0
            L1=fc22(1,2);
            L2=fc22(1,3);
            L3=fc22(1,4);
        else
            L1=fc22(1,3);
            L2=fc22(1,2);
            L3=fc22(1,1);
        end
        
        %% S{15}
        s15(k,1)=(L1+L2)/(L3-L2);
        
        %% S{16~17}
        s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
        s17(k,1)=-L2/(L3-L2);
        
        %% S{18~20}
        s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
        s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
        s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
    end
    S{15}=s15;
    S{16}=s16;
    S{17}=s17;
    S{18}=s18;
    S{19}=s19;
    S{20}=s20;
    
    %% S{21~26}
    S{21}=1./(S{5}.*S{6})';
    S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
    S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
    
    S{24}=(1./S{5}+1./S{6})';
    S{25}=S{6}'.*S{22};
    S{26}=S{5}'.*S{23};
    
    %% S{27~29}
    S{27}=(S{5}./(S{6}-S{5}))';
    S{28}=(S{4}.*S{7})';
    S{29}=S{8}';
    %% Fill
    FC{46}=S;
    FACE{l}=FC;
end

% FC{47}
for l=1:O
    FC=FACE{l};
    fc16=FC{16};
    fc22=FC{22};
    L_V_2=FC{4}*V2;
    
    %% initialization
    S=cell(1,20);
    s1=zeros(1,q2,'single');
    s2=zeros(1,q2,'single');
    s3=zeros(1,q2,'single');
    
    s4=zeros(1,q2);
    s5=zeros(1,q2);
    s6=zeros(1,q2);
    
    s7=zeros(1,q2);
    s8=zeros(1,q2);
    
    s9=zeros(1,q2);
    s10=zeros(1,q2);
    
    s11=zeros(1,q2);
    s12=zeros(1,q2);
    
    s13=zeros(1,q2);
    s14=zeros(1,q2);
    
    s15=zeros(1,q2)';
    
    s16=zeros(1,q2)';
    s17=zeros(1,q2)';
    
    s18=zeros(1,q2)';
    s19=zeros(1,q2)';
    s20=zeros(1,q2)';
    
    %% S{1}, S{2} and S{3}
    for k=1:q2
        if L_V_2(k)>=0
            s1(1,k)=fc16(1,2);
            s2(1,k)=fc16(1,3);
            s3(1,k)=fc16(1,4);
        else
            s1(1,k)=fc16(1,3);
            s2(1,k)=fc16(1,2);
            s3(1,k)=fc16(1,1);
        end
    end
    S{1}=s1;
    S{2}=s2;
    S{3}=s3;
    
    %% S{4}, S{5} and S{6}
    for k=1:q2
        if L_V_2(k)>=0
            s4(1,k)=fc22(1,2);
            s5(1,k)=fc22(1,3);
            s6(1,k)=fc22(1,4);
        else
            s4(1,k)=fc22(1,3);
            s5(1,k)=fc22(1,2);
            s6(1,k)=fc22(1,1);
        end
    end
    S{4}=s4;
    S{5}=s5;
    S{6}=s6;
    
    %% S{7~8}
    for k=1:q2
        s7(1,k)=fc22(1,2)+fc22(1,3);
        if L_V_2(k)>=0
            s8(1,k)=fc22(1,3)/s7(1,k);
        else
            s8(1,k)=fc22(1,2)/s7(1,k);
        end
    end
    S{7}=1./s7;
    S{8}=s8;
    
    %% S{9~10}
    for k=1:q2
        if L_V_2(k)>=0
            s9(1,k)=fc22(1,4)-fc22(1,3);
            s10(1,k)=fc22(1,3)/s9(1,k);
        else
            s9(1,k)=fc22(1,1)-fc22(1,2);
            s10(1,k)=fc22(1,2)/s9(1,k);
        end
    end
    S{9}=1./s9;
    S{10}=s10;
    
    %% S{11~12}
    for k=1:q2
        if L_V_2(k)>=0
            s11(1,k)=fc22(1,2)+fc22(1,4);
            s12(1,k)=fc22(1,3)/s11(1,k);
        else
            s11(1,k)=fc22(1,1)+fc22(1,3);
            s12(1,k)=fc22(1,2)/s11(1,k);
        end
    end
    S{11}=1./s11;
    S{12}=s12;
    
    %% S{13~14}
    for k=1:q2
        if L_V_2(k)>=0
            s13(1,k)=fc22(1,3);
        else
            s13(1,k)=fc22(1,2);
        end
    end
    S{13}=1./s13;
    S{14}=s14+1;
    
    %% S{15~20}
    for k=1:q2
        if L_V_2(k)>=0
            L1=fc22(1,2);
            L2=fc22(1,3);
            L3=fc22(1,4);
        else
            L1=fc22(1,3);
            L2=fc22(1,2);
            L3=fc22(1,1);
        end
        
        %% S{15}
        s15(k,1)=(L1+L2)/(L3-L2);
        
        %% S{16~17}
        s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
        s17(k,1)=-L2/(L3-L2);
        
        %% S{18~20}
        s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
        s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
        s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
    end
    S{15}=s15;
    S{16}=s16;
    S{17}=s17;
    S{18}=s18;
    S{19}=s19;
    S{20}=s20;
    
    %% S{21~26}
    S{21}=1./(S{5}.*S{6})';
    S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
    S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
    
    S{24}=(1./S{5}+1./S{6})';
    S{25}=S{6}'.*S{22};
    S{26}=S{5}'.*S{23};
    
    %% S{27~29}
    S{27}=(S{5}./(S{6}-S{5}))';
    S{28}=(S{4}.*S{7})';
    S{29}=S{8}';
    %% Fill
    FC{47}=S;
    FACE{l}=FC;
end


%% FC{48}, FC{49}
% FC{48}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc25=FC{25};
        fc27=FC{27};
        L_V_1=FC{4}*V1;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q1,'single');
        s2=zeros(1,q1,'single');
        s3=zeros(1,q1,'single');
        
        s4=zeros(1,q1);
        s5=zeros(1,q1);
        s6=zeros(1,q1);
        
        s7=zeros(1,q1);
        s8=zeros(1,q1);
        
        s9=zeros(1,q1);
        s10=zeros(1,q1);
        
        s11=zeros(1,q1);
        s12=zeros(1,q1);
        
        s13=zeros(1,q1);
        s14=zeros(1,q1);
        
        s15=zeros(1,q1)';
        
        s16=zeros(1,q1)';
        s17=zeros(1,q1)';
        
        s18=zeros(1,q1)';
        s19=zeros(1,q1)';
        s20=zeros(1,q1)';
        %% S{1}, S{2} and S{3}
        for k=1:q1
            if L_V_1(k)>=0
                s1(1,k)=fc25(1,2);
                s2(1,k)=fc25(1,3);
                s3(1,k)=fc25(1,4);
            else
                s1(1,k)=fc25(1,3);
                s2(1,k)=fc25(1,2);
                s3(1,k)=fc25(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q1
            if L_V_1(k)>=0
                s4(1,k)=fc27(1,2);
                s5(1,k)=fc27(1,3);
                s6(1,k)=fc27(1,4);
            else
                s4(1,k)=fc27(1,3);
                s5(1,k)=fc27(1,2);
                s6(1,k)=fc27(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q1
            s7(1,k)=fc27(1,2)+fc27(1,3);
            if L_V_1(k)>=0
                s8(1,k)=fc27(1,3)/s7(1,k);
            else
                s8(1,k)=fc27(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q1
            if L_V_1(k)>=0
                s9(1,k)=fc27(1,4)-fc27(1,3);
                s10(1,k)=fc27(1,3)/s9(1,k);
            else
                s9(1,k)=fc27(1,1)-fc27(1,2);
                s10(1,k)=fc27(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q1
            if L_V_1(k)>=0
                s11(1,k)=fc27(1,2)+fc27(1,4);
                s12(1,k)=fc27(1,3)/s11(1,k);
            else
                s11(1,k)=fc27(1,1)+fc27(1,3);
                s12(1,k)=fc27(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q1
            if L_V_1(k)>=0
                s13(1,k)=fc27(1,3);
            else
                s13(1,k)=fc27(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q1
            if L_V_1(k)>=0
                L1=fc27(1,2);
                L2=fc27(1,3);
                L3=fc27(1,4);
            else
                L1=fc27(1,3);
                L2=fc27(1,2);
                L3=fc27(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{48}=S;
        FACE{l}=FC;
    end
end

% FC{49}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc25=FC{25};
        fc27=FC{27};
        L_V_2=FC{4}*V2;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q2,'single');
        s2=zeros(1,q2,'single');
        s3=zeros(1,q2,'single');
        
        s4=zeros(1,q2);
        s5=zeros(1,q2);
        s6=zeros(1,q2);
        
        s7=zeros(1,q2);
        s8=zeros(1,q2);
        
        s9=zeros(1,q2);
        s10=zeros(1,q2);
        
        s11=zeros(1,q2);
        s12=zeros(1,q2);
        
        s13=zeros(1,q2);
        s14=zeros(1,q2);
        
        s15=zeros(1,q2)';
        
        s16=zeros(1,q2)';
        s17=zeros(1,q2)';
        
        s18=zeros(1,q2)';
        s19=zeros(1,q2)';
        s20=zeros(1,q2)';
        %% S{1}, S{2} and S{3}
        for k=1:q2
            if L_V_2(k)>=0
                s1(1,k)=fc25(1,2);
                s2(1,k)=fc25(1,3);
                s3(1,k)=fc25(1,4);
            else
                s1(1,k)=fc25(1,3);
                s2(1,k)=fc25(1,2);
                s3(1,k)=fc25(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q2
            if L_V_2(k)>=0
                s4(1,k)=fc27(1,2);
                s5(1,k)=fc27(1,3);
                s6(1,k)=fc27(1,4);
            else
                s4(1,k)=fc27(1,3);
                s5(1,k)=fc27(1,2);
                s6(1,k)=fc27(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q2
            s7(1,k)=fc27(1,2)+fc27(1,3);
            if L_V_2(k)>=0
                s8(1,k)=fc27(1,3)/s7(1,k);
            else
                s8(1,k)=fc27(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q2
            if L_V_2(k)>=0
                s9(1,k)=fc27(1,4)-fc27(1,3);
                s10(1,k)=fc27(1,3)/s9(1,k);
            else
                s9(1,k)=fc27(1,1)-fc27(1,2);
                s10(1,k)=fc27(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q2
            if L_V_2(k)>=0
                s11(1,k)=fc27(1,2)+fc27(1,4);
                s12(1,k)=fc27(1,3)/s11(1,k);
            else
                s11(1,k)=fc27(1,1)+fc27(1,3);
                s12(1,k)=fc27(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q2
            if L_V_2(k)>=0
                s13(1,k)=fc27(1,3);
            else
                s13(1,k)=fc27(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q2
            if L_V_2(k)>=0
                L1=fc27(1,2);
                L2=fc27(1,3);
                L3=fc27(1,4);
            else
                L1=fc27(1,3);
                L2=fc27(1,2);
                L3=fc27(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{49}=S;
        FACE{l}=FC;
    end
end


%% FC{50}, FC{51}
% FC{50}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc31=FC{31};
        fc33=FC{33};
        L_V_1=FC{4}*V1;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q1,'single');
        s2=zeros(1,q1,'single');
        s3=zeros(1,q1,'single');
        
        s4=zeros(1,q1);
        s5=zeros(1,q1);
        s6=zeros(1,q1);
        
        s7=zeros(1,q1);
        s8=zeros(1,q1);
        
        s9=zeros(1,q1);
        s10=zeros(1,q1);
        
        s11=zeros(1,q1);
        s12=zeros(1,q1);
        
        s13=zeros(1,q1);
        s14=zeros(1,q1);
        
        s15=zeros(1,q1)';
        
        s16=zeros(1,q1)';
        s17=zeros(1,q1)';
        
        s18=zeros(1,q1)';
        s19=zeros(1,q1)';
        s20=zeros(1,q1)';
        %% S{1}, S{2} and S{3}
        for k=1:q1
            if L_V_1(k)>=0
                s1(1,k)=fc31(1,2);
                s2(1,k)=fc31(1,3);
                s3(1,k)=fc31(1,4);
            else
                s1(1,k)=fc31(1,3);
                s2(1,k)=fc31(1,2);
                s3(1,k)=fc31(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q1
            if L_V_1(k)>=0
                s4(1,k)=fc33(1,2);
                s5(1,k)=fc33(1,3);
                s6(1,k)=fc33(1,4);
            else
                s4(1,k)=fc33(1,3);
                s5(1,k)=fc33(1,2);
                s6(1,k)=fc33(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q1
            s7(1,k)=fc33(1,2)+fc33(1,3);
            if L_V_1(k)>=0
                s8(1,k)=fc33(1,3)/s7(1,k);
            else
                s8(1,k)=fc33(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q1
            if L_V_1(k)>=0
                s9(1,k)=fc33(1,4)-fc33(1,3);
                s10(1,k)=fc33(1,3)/s9(1,k);
            else
                s9(1,k)=fc33(1,1)-fc33(1,2);
                s10(1,k)=fc33(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q1
            if L_V_1(k)>=0
                s11(1,k)=fc33(1,2)+fc33(1,4);
                s12(1,k)=fc33(1,3)/s11(1,k);
            else
                s11(1,k)=fc33(1,1)+fc33(1,3);
                s12(1,k)=fc33(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q1
            if L_V_1(k)>=0
                s13(1,k)=fc33(1,3);
            else
                s13(1,k)=fc33(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q1
            if L_V_1(k)>=0
                L1=fc33(1,2);
                L2=fc33(1,3);
                L3=fc33(1,4);
            else
                L1=fc33(1,3);
                L2=fc33(1,2);
                L3=fc33(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{50}=S;
        FACE{l}=FC;
    end
end

% FC{51}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc31=FC{31};
        fc33=FC{33};
        L_V_2=FC{4}*V2;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q2,'single');
        s2=zeros(1,q2,'single');
        s3=zeros(1,q2,'single');
        
        s4=zeros(1,q2);
        s5=zeros(1,q2);
        s6=zeros(1,q2);
        
        s7=zeros(1,q2);
        s8=zeros(1,q2);
        
        s9=zeros(1,q2);
        s10=zeros(1,q2);
        
        s11=zeros(1,q2);
        s12=zeros(1,q2);
        
        s13=zeros(1,q2);
        s14=zeros(1,q2);
        
        s15=zeros(1,q2)';
        
        s16=zeros(1,q2)';
        s17=zeros(1,q2)';
        
        s18=zeros(1,q2)';
        s19=zeros(1,q2)';
        s20=zeros(1,q2)';
        %% S{1}, S{2} and S{3}
        for k=1:q2
            if L_V_2(k)>=0
                s1(1,k)=fc31(1,2);
                s2(1,k)=fc31(1,3);
                s3(1,k)=fc31(1,4);
            else
                s1(1,k)=fc31(1,3);
                s2(1,k)=fc31(1,2);
                s3(1,k)=fc31(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q2
            if L_V_2(k)>=0
                s4(1,k)=fc33(1,2);
                s5(1,k)=fc33(1,3);
                s6(1,k)=fc33(1,4);
            else
                s4(1,k)=fc33(1,3);
                s5(1,k)=fc33(1,2);
                s6(1,k)=fc33(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q2
            s7(1,k)=fc33(1,2)+fc33(1,3);
            if L_V_2(k)>=0
                s8(1,k)=fc33(1,3)/s7(1,k);
            else
                s8(1,k)=fc33(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q2
            if L_V_2(k)>=0
                s9(1,k)=fc33(1,4)-fc33(1,3);
                s10(1,k)=fc33(1,3)/s9(1,k);
            else
                s9(1,k)=fc33(1,1)-fc33(1,2);
                s10(1,k)=fc33(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q2
            if L_V_2(k)>=0
                s11(1,k)=fc33(1,2)+fc33(1,4);
                s12(1,k)=fc33(1,3)/s11(1,k);
            else
                s11(1,k)=fc33(1,1)+fc33(1,3);
                s12(1,k)=fc33(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q2
            if L_V_2(k)>=0
                s13(1,k)=fc33(1,3);
            else
                s13(1,k)=fc33(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q2
            if L_V_2(k)>=0
                L1=fc33(1,2);
                L2=fc33(1,3);
                L3=fc33(1,4);
            else
                L1=fc33(1,3);
                L2=fc33(1,2);
                L3=fc33(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{51}=S;
        FACE{l}=FC;
    end
end


%% FC{52}, FC{53}
% FC{52}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc38=FC{38};
        fc40=FC{40};
        L_V_1=FC{4}*V1;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q1,'single');
        s2=zeros(1,q1,'single');
        s3=zeros(1,q1,'single');
        
        s4=zeros(1,q1);
        s5=zeros(1,q1);
        s6=zeros(1,q1);
        
        s7=zeros(1,q1);
        s8=zeros(1,q1);
        
        s9=zeros(1,q1);
        s10=zeros(1,q1);
        
        s11=zeros(1,q1);
        s12=zeros(1,q1);
        
        s13=zeros(1,q1);
        s14=zeros(1,q1);
        
        s15=zeros(1,q1)';
        
        s16=zeros(1,q1)';
        s17=zeros(1,q1)';
        
        s18=zeros(1,q1)';
        s19=zeros(1,q1)';
        s20=zeros(1,q1)';
        %% S{1}, S{2} and S{3}
        for k=1:q1
            if L_V_1(k)>=0
                s1(1,k)=fc38(1,2);
                s2(1,k)=fc38(1,3);
                s3(1,k)=fc38(1,4);
            else
                s1(1,k)=fc38(1,3);
                s2(1,k)=fc38(1,2);
                s3(1,k)=fc38(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q1
            if L_V_1(k)>=0
                s4(1,k)=fc40(1,2);
                s5(1,k)=fc40(1,3);
                s6(1,k)=fc40(1,4);
            else
                s4(1,k)=fc40(1,3);
                s5(1,k)=fc40(1,2);
                s6(1,k)=fc40(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q1
            s7(1,k)=fc40(1,2)+fc40(1,3);
            if L_V_1(k)>=0
                s8(1,k)=fc40(1,3)/s7(1,k);
            else
                s8(1,k)=fc40(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q1
            if L_V_1(k)>=0
                s9(1,k)=fc40(1,4)-fc40(1,3);
                s10(1,k)=fc40(1,3)/s9(1,k);
            else
                s9(1,k)=fc40(1,1)-fc40(1,2);
                s10(1,k)=fc40(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q1
            if L_V_1(k)>=0
                s11(1,k)=fc40(1,2)+fc40(1,4);
                s12(1,k)=fc40(1,3)/s11(1,k);
            else
                s11(1,k)=fc40(1,1)+fc40(1,3);
                s12(1,k)=fc40(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q1
            if L_V_1(k)>=0
                s13(1,k)=fc40(1,3);
            else
                s13(1,k)=fc40(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q1
            if L_V_1(k)>=0
                L1=fc40(1,2);
                L2=fc40(1,3);
                L3=fc40(1,4);
            else
                L1=fc40(1,3);
                L2=fc40(1,2);
                L3=fc40(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{52}=S;
        FACE{l}=FC;
    end
end

% FC{53}
for l=1:O
    FC=FACE{l};
    if FC{24}==2 || FC{24}==3
        fc38=FC{38};
        fc40=FC{40};
        L_V_2=FC{4}*V2;
        
        %% initialization
        S=cell(1,20);
        s1=zeros(1,q2,'single');
        s2=zeros(1,q2,'single');
        s3=zeros(1,q2,'single');
        
        s4=zeros(1,q2);
        s5=zeros(1,q2);
        s6=zeros(1,q2);
        
        s7=zeros(1,q2);
        s8=zeros(1,q2);
        
        s9=zeros(1,q2);
        s10=zeros(1,q2);
        
        s11=zeros(1,q2);
        s12=zeros(1,q2);
        
        s13=zeros(1,q2);
        s14=zeros(1,q2);
        
        s15=zeros(1,q2)';
        
        s16=zeros(1,q2)';
        s17=zeros(1,q2)';
        
        s18=zeros(1,q2)';
        s19=zeros(1,q2)';
        s20=zeros(1,q2)';
        %% S{1}, S{2} and S{3}
        for k=1:q2
            if L_V_2(k)>=0
                s1(1,k)=fc38(1,2);
                s2(1,k)=fc38(1,3);
                s3(1,k)=fc38(1,4);
            else
                s1(1,k)=fc38(1,3);
                s2(1,k)=fc38(1,2);
                s3(1,k)=fc38(1,1);
            end
        end
        S{1}=s1;
        S{2}=s2;
        S{3}=s3;
        
        %% S{4}, S{5} and S{6}
        for k=1:q2
            if L_V_2(k)>=0
                s4(1,k)=fc40(1,2);
                s5(1,k)=fc40(1,3);
                s6(1,k)=fc40(1,4);
            else
                s4(1,k)=fc40(1,3);
                s5(1,k)=fc40(1,2);
                s6(1,k)=fc40(1,1);
            end
        end
        S{4}=s4;
        S{5}=s5;
        S{6}=s6;
        
        %% S{7~8}
        for k=1:q2
            s7(1,k)=fc40(1,2)+fc40(1,3);
            if L_V_2(k)>=0
                s8(1,k)=fc40(1,3)/s7(1,k);
            else
                s8(1,k)=fc40(1,2)/s7(1,k);
            end
        end
        S{7}=1./s7;
        S{8}=s8;
        
        %% S{9~10}
        for k=1:q2
            if L_V_2(k)>=0
                s9(1,k)=fc40(1,4)-fc40(1,3);
                s10(1,k)=fc40(1,3)/s9(1,k);
            else
                s9(1,k)=fc40(1,1)-fc40(1,2);
                s10(1,k)=fc40(1,2)/s9(1,k);
            end
        end
        S{9}=1./s9;
        S{10}=s10;
        
        %% S{11~12}
        for k=1:q2
            if L_V_2(k)>=0
                s11(1,k)=fc40(1,2)+fc40(1,4);
                s12(1,k)=fc40(1,3)/s11(1,k);
            else
                s11(1,k)=fc40(1,1)+fc40(1,3);
                s12(1,k)=fc40(1,2)/s11(1,k);
            end
        end
        S{11}=1./s11;
        S{12}=s12;
        
        %% S{13~14}
        for k=1:q2
            if L_V_2(k)>=0
                s13(1,k)=fc40(1,3);
            else
                s13(1,k)=fc40(1,2);
            end
        end
        S{13}=1./s13;
        S{14}=s14+1;
        
        %% S{15~20}
        for k=1:q2
            if L_V_2(k)>=0
                L1=fc40(1,2);
                L2=fc40(1,3);
                L3=fc40(1,4);
            else
                L1=fc40(1,3);
                L2=fc40(1,2);
                L3=fc40(1,1);
            end
            
            %% S{15}
            s15(k,1)=(L1+L2)/(L3-L2);
            
            %% S{16~17}
            s16(k,1)=(L1+L2*(L1+L3)/(L3-L2))/(L1+L2);
            s17(k,1)=-L2/(L3-L2);
            
            %% S{18~20}
            s18(k,1)=L2*L3/(L1+L2)/(L1+L3);
            s19(k,1)=(1+L2/(L3-L2))*L1/(L1+L2);
            s20(k,1)=-L1*L2/(L1+L3)/(L3-L2);
        end
        S{15}=s15;
        S{16}=s16;
        S{17}=s17;
        S{18}=s18;
        S{19}=s19;
        S{20}=s20;
        
        %% S{21~26}
        S{21}=1./(S{5}.*S{6})';
        S{22}=1./(S{5}.*S{5}-S{5}.*S{6})';
        S{23}=1./(S{6}.*S{6}-S{5}.*S{6})';
        
        S{24}=(1./S{5}+1./S{6})';
        S{25}=S{6}'.*S{22};
        S{26}=S{5}'.*S{23};
        
        %% S{27~29}
        S{27}=(S{5}./(S{6}-S{5}))';
        S{28}=(S{4}.*S{7})';
        S{29}=S{8}';
        %% Fill
        FC{53}=S;
        FACE{l}=FC;
    end
end


%% FC{54}, FC{55}
% FC{54}
for l=1:O
    FC=FACE{l};
    L_V_1=FC{4}*V1;
    fc54=zeros(q1,1);
    for k=1:q1
        if L_V_1(k)>=0
            fc54(k,1)=1;
        else
            fc54(k,1)=0;
        end
    end
    %% Fill
    FC{54}=fc54;
    FACE{l}=FC;
end

% FC{55}
for l=1:O
    FC=FACE{l};
    L_V_2=FC{4}*V2;
    fc55=zeros(q2,1);
    for k=1:q2
        if L_V_2(k)>=0
            fc55(k,1)=1;
        else
            fc55(k,1)=0;
        end
    end
    %% Fill
    FC{55}=fc55;
    FACE{l}=FC;
end


%% The following section is the old code for filling the data in FACE before v.11.24.2015      
    %     %% Filling FC{16}, {17}
    %     for l=1:O
    %         FC=FACE{l};
    %         S_V1=cell(K,1); % stencil for V1
    %         S_V2=cell(K,1); % stencil for V2
    %         neigh_up=FC{12};
    %         neigh_down=FC{13};
    %         %% S{1} and S{2}
    %         % S{1}, S{2} for V1
    %         DC1=zeros(1,q1);
    %         UC1=zeros(1,q1);
    %         for k=1:q1
    %             if FC{4}*V1(:,k)>=0
    %                 DC1(1,k)=neigh_down(1,1);
    %                 UC1(1,k)=neigh_up(1,1);
    %             else
    %                 DC1(1,k)=neigh_up(1,1);
    %                 UC1(1,k)=neigh_down(1,1);
    %             end
    %         end
    %         % Check
    %         % Cell pair
    %         for k=1:q1
    %             if length(unique(union([DC1(1,k),UC1(1,k)],[neigh_down(1,1),neigh_up(1,1)])))~=length([neigh_down(1,1),neigh_up(1,1)])
    %                 error('The downwind & upwind cells are not correct!')
    %             end
    %         end
    %         % Direction
    %         for k=1:q1
    %             if DC1(1,k)~=0 && UC1(1,k)~=0
    %                 Cell_down=CELL{DC1(1,k)};
    %                 Cell_up=CELL{UC1(1,k)};
    %                 n_u2d=(Cell_down{5}-Cell_up{5})';
    %                 if single(e+n_u2d*V1(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             elseif DC1(1,k)==0 && UC1(1,k)~=0
    %                 Cell_up=CELL{UC1(1,k)};
    %                 n_u2d=(FC{7}-Cell_up{5})';
    %                 if single(e+n_u2d*V1(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             elseif DC1(1,k)~=0 && UC1(1,k)==0
    %                 Cell_down=CELL{DC1(1,k)};
    %                 n_u2d=(Cell_down{5}-FC{7})';
    %                 if single(e+n_u2d*V1(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             else
    %                 error('There is error in FC{12} and FC{13}!');
    %             end
    %         end
    %         % S{1}, S{2} for V2
    %         DC2=zeros(1,q2);
    %         UC2=zeros(1,q2);
    %         for k=1:q2
    %             if FC{4}*V2(:,k)>=0
    %                 DC2(1,k)=neigh_down(1,1);
    %                 UC2(1,k)=neigh_up(1,1);
    %             else
    %                 DC2(1,k)=neigh_up(1,1);
    %                 UC2(1,k)=neigh_down(1,1);
    %             end
    %         end
    %         % Check
    %         % Cell pair
    %         for k=1:q2
    %             if length(unique(union([DC2(1,k),UC2(1,k)],[neigh_down(1,1),neigh_up(1,1)])))~=length([neigh_down(1,1),neigh_up(1,1)])
    %                 error('The downwind & upwind cells are not correct!')
    %             end
    %         end
    %         % Direction
    %         for k=1:q2
    %             if DC2(1,k)~=0 && UC2(1,k)~=0
    %                 Cell_down=CELL{DC2(1,k)};
    %                 Cell_up=CELL{UC2(1,k)};
    %                 n_u2d=(Cell_down{5}-Cell_up{5})';
    %                 if single(e+n_u2d*V2(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             elseif DC2(1,k)==0 && UC2(1,k)~=0
    %                 Cell_up=CELL{UC2(1,k)};
    %                 n_u2d=(FC{7}-Cell_up{5})';
    %                 if single(e+n_u2d*V2(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             elseif DC2(1,k)~=0 && UC2(1,k)==0
    %                 Cell_down=CELL{DC2(1,k)};
    %                 n_u2d=(Cell_down{5}-FC{7})';
    %                 if single(e+n_u2d*V2(:,k))<single(e)
    %                     error('The downwind & upwind cells should be switched!');
    %                 end
    %             else
    %                 error('There is error in FC{12} and FC{13}!');
    %             end
    %         end
    %         % Fill the data
    %         S_V1{1,1}=DC1;
    %         S_V1{2,1}=UC1;
    %
    %         S_V2{1,1}=DC2;
    %         S_V2{2,1}=UC2;
    %         %% S{4}~S{7}
    %         %%%%%%% The output of this section C_b, C_j, C_j_down, C_j_up, d_cc, d_nc, dc_down, dc_up can be used multiple times
    %         if FC{2}~=0 % current face is on boundary
    %             C_b=FC{7};
    %             neigh_up=FC{12};
    %             neigh_down=FC{13};
    %             cell_in=setxor(0,[neigh_down(1,1),neigh_up(1,1)]);
    %             if length(cell_in)~=1
    %                 error('There is only one cell attached to the boundary face!');
    %             end
    %             P=CELL{cell_in};
    %             C_j=norm_joint(C_b,FC{4}',P{5});
    %             if in_triangle(C_j,P{13},P{14},P{15})
    %                 ;
    %             else
    %                 error('The generated point is not within the current cell!');
    %             end
    %             d_cc=dis(C_b,C_j); % distance from the face midpoint
    %             % to the point that is closest to the centroid of current cell
    %             d_nc=d_cc; % distance from the face midpoint
    %             % to the point that is closest to the centroid of neighbor cell
    %         else  % Internal edge
    %             C_b=FC{7};
    %             neigh_up=FC{12};
    %             neigh_down=FC{13};
    %             P=CELL{neigh_down(1,1)};
    %             Q=CELL{neigh_up(1,1)};
    %             C_j_down=norm_joint(C_b,FC{4}',P{5});
    %             C_j_up=norm_joint(C_b,FC{4}',Q{5});
    %             if in_triangle(C_j_down,P{13},P{14},P{15})
    %                 ;
    %             else
    %                 error('The generated point is not within the current cell!');
    %             end
    %             if in_triangle(C_j_up,Q{13},Q{14},Q{15})
    %                 ;
    %             else
    %                 error('The generated point is not within the current cell!');
    %             end
    %             d_c_down=dis(C_b,C_j_down);
    %             d_c_up=dis(C_b,C_j_up);
    %         end
    %         %%%%%%% The output of this section C_b, C_j, C_j_down, C_j_up, d_cc, d_nc, dc_down, dc_up can be used multiple times
    %         % S{4}~S{7} for V1
    %         DC1=S_V1{1,1};
    %         UC1=S_V1{2,1};
    %         x_down1=zeros(1,q1);
    %         y_down1=zeros(1,q1);
    %         x_up1=zeros(1,q1);
    %         y_up1=zeros(1,q1);
    %         for k=1:q1
    %             if DC1(1,k)==neigh_down(1,1) && UC1(1,k)==neigh_up(1,1)
    %                 if FC{2}~=0
    %                     if DC1(1,k)==0 && UC1(1,k)~=0
    %                         x_up1(1,k)=C_j(1,1);
    %                         y_up1(1,k)=C_j(2,1);
    %                         x_down1(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_down1(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_down1(1,k)<X2 && x_down1(1,k)>X1) && (y_down1(1,k)<Y2 && y_down1(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     elseif DC1(1,k)~=0 && UC1(1,k)==0
    %                         x_down1(1,k)=C_j(1,1);
    %                         y_down1(1,k)=C_j(2,1);
    %                         x_up1(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_up1(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_up1(1,k)<X2 && x_up1(1,k)>X1) && (y_up1(1,k)<Y2 && y_up1(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     else
    %                         error('The face on the boundary contains false neighbor cell info!');
    %                     end
    %                 else
    %                     x_down1(1,k)=C_j_down(1,1);
    %                     y_down1(1,k)=C_j_down(2,1);
    %                     x_up1(1,k)=C_j_up(1,1);
    %                     y_up1(1,k)=C_j_up(2,1);
    %                     if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                     if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                 end
    %             elseif DC1(1,k)==neigh_up(1,1) && UC1(1,k)==neigh_down(1,1)
    %                 if FC{2}~=0
    %                     if DC1(1,k)==0 && UC1(1,k)~=0
    %                         x_up1(1,k)=C_j(1,1);
    %                         y_up1(1,k)=C_j(2,1);
    %                         x_down1(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_down1(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_down1(1,k)<X2 && x_down1(1,k)>X1) && (y_down1(1,k)<Y2 && y_down1(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     elseif DC1(1,k)~=0 && UC1(1,k)==0
    %                         x_down1(1,k)=C_j(1,1);
    %                         y_down1(1,k)=C_j(2,1);
    %                         x_up1(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_up1(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_up1(1,k)<X2 && x_up1(1,k)>X1) && (y_up1(1,k)<Y2 && y_up1(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     else
    %                         error('The face on the boundary contains false neighbor cell info!');
    %                     end
    %                 else
    %                     x_down1(1,k)=C_j_up(1,1);
    %                     y_down1(1,k)=C_j_up(2,1);
    %                     x_up1(1,k)=C_j_down(1,1);
    %                     y_up1(1,k)=C_j_down(2,1);
    %                     if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                     if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                 end
    %             else
    %                 error('Two neighbor cells for the current face reside both on one side of the face!');
    %             end
    %         end
    %         % Fill data S{4}~S{7} for V1
    %         S_V1{4,1}=x_down1;
    %         S_V1{5,1}=y_down1;
    %         S_V1{6,1}=x_up1;
    %         S_V1{7,1}=y_up1;
    %
    %         % S{4}~S{7} for V2
    %         DC2=S_V2{1,1};
    %         UC2=S_V2{2,1};
    %         x_down2=zeros(1,q2);
    %         y_down2=zeros(1,q2);
    %         x_up2=zeros(1,q2);
    %         y_up2=zeros(1,q2);
    %         for k=1:q2
    %             if DC2(1,k)==neigh_down(1,1) && UC2(1,k)==neigh_up(1,1)
    %                 if FC{2}~=0
    %                     if DC2(1,k)==0 && UC2(1,k)~=0
    %                         x_up2(1,k)=C_j(1,1);
    %                         y_up2(1,k)=C_j(2,1);
    %                         x_down2(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_down2(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_down2(1,k)<X2 && x_down2(1,k)>X1) && (y_down2(1,k)<Y2 && y_down2(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     elseif DC2(1,k)~=0 && UC2(1,k)==0
    %                         x_down2(1,k)=C_j(1,1);
    %                         y_down2(1,k)=C_j(2,1);
    %                         x_up2(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_up2(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_up2(1,k)<X2 && x_up2(1,k)>X1) && (y_up2(1,k)<Y2 && y_up2(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     else
    %                         error('The face on the boundary contains false neighbor cell info!');
    %                     end
    %                 else
    %                     x_down2(1,k)=C_j_down(1,1);
    %                     y_down2(1,k)=C_j_down(2,1);
    %                     x_up2(1,k)=C_j_up(1,1);
    %                     y_up2(1,k)=C_j_up(2,1);
    %                     if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                     if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                 end
    %             elseif DC2(1,k)==neigh_up(1,1) && UC2(1,k)==neigh_down(1,1)
    %                 if FC{2}~=0
    %                     if DC2(1,k)==0 && UC2(1,k)~=0
    %                         x_up2(1,k)=C_j(1,1);
    %                         y_up2(1,k)=C_j(2,1);
    %                         x_down2(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_down2(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_down2(1,k)<X2 && x_down2(1,k)>X1) && (y_down2(1,k)<Y2 && y_down2(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     elseif DC2(1,k)~=0 && UC2(1,k)==0
    %                         x_down2(1,k)=C_j(1,1);
    %                         y_down2(1,k)=C_j(2,1);
    %                         x_up2(1,k)=2*C_b(1,1)-C_j(1,1);
    %                         y_up2(1,k)=2*C_b(2,1)-C_j(2,1);
    %                         if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
    %                             error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
    %                         end
    %                         if (x_up2(1,k)<X2 && x_up2(1,k)>X1) && (y_up2(1,k)<Y2 && y_up2(1,k)>Y1)
    %                             error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
    %                         end
    %                     else
    %                         error('The face on the boundary contains false neighbor cell info!');
    %                     end
    %                 else
    %                     x_down2(1,k)=C_j_up(1,1);
    %                     y_down2(1,k)=C_j_up(2,1);
    %                     x_up2(1,k)=C_j_down(1,1);
    %                     y_up2(1,k)=C_j_down(2,1);
    %                     if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                     if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
    %                         error('The stencil point for any side of the interior face should be located inside of the computational domain!');
    %                     end
    %                 end
    %             else
    %                 error('Two neighbor cells for the current face reside both on one side of the face!');
    %             end
    %         end
    %         % Fill data S{4}~S{7} for V2
    %         S_V2{4,1}=x_down2;
    %         S_V2{5,1}=y_down2;
    %         S_V2{6,1}=x_up2;
    %         S_V2{7,1}=y_up2;
    %
    %         %% S{13}, S{15}
    %         % S{13}, S{15} for V1
    %         x_down1=S_V1{4,1};
    %         y_down1=S_V1{5,1};
    %         x_up1=S_V1{6,1};
    %         y_up1=S_V1{7,1};
    %         dis_u2d=zeros(1,q1);
    %         One_over_dis_u2d=zeros(1,q1);
    %         for k=1:q1
    %             dis_u2d(1,k)=dis([x_down1(1,k);y_down1(1,k)],[x_up1(1,k);y_up1(1,k)]);
    %             if FC{2}~=0
    %                 if single(e+dis_u2d(1,k))~=single(e+2*d_nc)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             else
    %                 if single(e+dis_u2d(1,k))~=single(e+d_c_down+d_c_up)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             end
    %             One_over_dis_u2d(1,k)=1/dis_u2d(1,k);
    %         end
    %         % Fill the data S{13}, S{15} for V1
    %         S_V1{13,1}=dis_u2d;
    %         S_V1{15,1}=One_over_dis_u2d;
    %
    %         % S{13}, S{15} for V2
    %         x_down2=S_V2{4,1};
    %         y_down2=S_V2{5,1};
    %         x_up2=S_V2{6,1};
    %         y_up2=S_V2{7,1};
    %         dis_u2d=zeros(1,q2);
    %         One_over_dis_u2d=zeros(1,q2);
    %         for k=1:q2
    %             dis_u2d(1,k)=dis([x_down2(1,k);y_down2(1,k)],[x_up2(1,k);y_up2(1,k)]);
    %             if FC{2}~=0
    %                 if single(e+dis_u2d(1,k))~=single(e+2*d_nc)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             else
    %                 if single(e+dis_u2d(1,k))~=single(e+d_c_down+d_c_up)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             end
    %             One_over_dis_u2d(1,k)=1/dis_u2d(1,k);
    %         end
    %         % Fill the data S{13}, S{15} for V2
    %         S_V2{13,1}=dis_u2d;
    %         S_V2{15,1}=One_over_dis_u2d;
    %
    %         %% S{17}, S{19}
    %         % S{17}, S{19} for V1
    %         x_down1=S_V1{4,1};
    %         y_down1=S_V1{5,1};
    %         x_up1=S_V1{6,1};
    %         y_up1=S_V1{7,1};
    %         u2d_over_u=zeros(1,q1);
    %         One_u2d_over_u=zeros(1,q1);
    %         for k=1:q1
    %             u2d_over_u(1,k)=dis([x_down1(1,k);y_down1(1,k)],[x_up1(1,k);y_up1(1,k)])/dis(C_b,[x_up1(1,k);y_up1(1,k)]);
    %             if FC{2}~=0
    %                 if single(e+u2d_over_u(1,k))~=single(e+2)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             else
    %                 if single(e+u2d_over_u(1,k))~=single(e+(d_c_down+d_c_up)/dis(C_b,[x_up1(1,k);y_up1(1,k)]))
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             end
    %             One_u2d_over_u(1,k)=1/u2d_over_u(1,k);
    %         end
    %         % Fill the data S{13}, S{15} for V1
    %         S_V1{17,1}=u2d_over_u;
    %         S_V1{19,1}=One_u2d_over_u;
    %
    %         % S{17}, S{19} for V2
    %         x_down2=S_V2{4,1};
    %         y_down2=S_V2{5,1};
    %         x_up2=S_V2{6,1};
    %         y_up2=S_V2{7,1};
    %         u2d_over_u=zeros(1,q2);
    %         One_u2d_over_u=zeros(1,q2);
    %         for k=1:q2
    %             u2d_over_u(1,k)=dis([x_down2(1,k);y_down2(1,k)],[x_up2(1,k);y_up2(1,k)])/dis(C_b,[x_up2(1,k);y_up2(1,k)]);
    %             if FC{2}~=0
    %                 if single(e+u2d_over_u(1,k))~=single(e+2)
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             else
    %                 if single(e+u2d_over_u(1,k))~=single(e+(d_c_down+d_c_up)/dis(C_b,[x_up2(1,k);y_up2(1,k)]))
    %                     error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
    %                 end
    %             end
    %             One_u2d_over_u(1,k)=1/u2d_over_u(1,k);
    %         end
    %         % Fill the data S{13}, S{15} for V1
    %         S_V2{17,1}=u2d_over_u;
    %         S_V2{19,1}=One_u2d_over_u;
    %         %% FC{18}
    %         fc18=zeros(2,3);
    %         if FC{2}==0
    %             fc18(:,1)=C_j_down;
    %             fc18(:,2)=C_j_up;
    %         else
    %             if neigh_up(1,1)==0 && neigh_down(1,1)~=0
    %                 fc18(:,1)=C_j;
    %             elseif neigh_up(1,1)~=0 && neigh_down(1,1)==0
    %                 fc18(:,2)=C_j;
    %             else
    %                 error('There is error in FC{12} and FC{13}!');
    %             end
    %         end
    %         FC{18}=fc18;
    %         %% S{10}, S{11}
    %         neigh_up=FC{12};
    %         neigh_down=FC{13};
    %         ZONE_UP=-1;
    %         ZONE_DOWN=-1;
    %         if FC{2}~=0
    %             if neigh_down(1,1)==0 && neigh_up(1,1)~=0
    %                 Cell_up=CELL{neigh_up(1,1)};
    %                 nd1_up=Cell_up{13};
    %                 nd2_up=Cell_up{14};
    %                 nd3_up=Cell_up{15};
    %                 in_zone_one=in_triangle(C_j,Cell_up{5},nd1_up,nd2_up);
    %                 in_zone_two=in_triangle(C_j,Cell_up{5},nd2_up,nd3_up);
    %                 in_zone_three=in_triangle(C_j,Cell_up{5},nd3_up,nd1_up);
    %                 if in_zone_one
    %                     ZONE_UP=1; % Located in zone 1
    %                 elseif in_zone_two
    %                     ZONE_UP=2; % Located in zone 2
    %                 elseif in_zone_three
    %                     ZONE_UP=3; % Located in zone 3
    %                 else
    %                     on_edge_one=on_edge(C_j,Cell_up{5},nd1_up);
    %                     on_edge_two=on_edge(C_j,Cell_up{5},nd2_up);
    %                     on_edge_three=on_edge(C_j,Cell_up{5},nd3_up);
    %                     if on_edge_one && on_edge_two && on_edge_three
    %                         ZONE_UP=4;
    %                     elseif on_edge_one
    %                         ZONE_UP=1;
    %                     elseif on_edge_two
    %                         ZONE_UP=2;
    %                     elseif on_edge_three
    %                         ZONE_UP=3;
    %                     else
    %                         error('The zone for the stencil point is not found!');
    %                     end
    %                 end
    %                 ZONE_DOWN=4;
    %             elseif neigh_down(1,1)~=0 && neigh_up(1,1)==0
    %                 Cell_down=CELL{neigh_down(1,1)};
    %                 nd1_down=Cell_down{13};
    %                 nd2_down=Cell_down{14};
    %                 nd3_down=Cell_down{15};
    %                 in_zone_one=in_triangle(C_j,Cell_down{5},nd1_down,nd2_down);
    %                 in_zone_two=in_triangle(C_j,Cell_down{5},nd2_down,nd3_down);
    %                 in_zone_three=in_triangle(C_j,Cell_down{5},nd3_down,nd1_down);
    %                 if in_zone_one
    %                     ZONE_DOWN=1; % Located in zone 1
    %                 elseif in_zone_two
    %                     ZONE_DOWN=2; % Located in zone 2
    %                 elseif in_zone_three
    %                     ZONE_DOWN=3; % Located in zone 3
    %                 else
    %                     on_edge_one=on_edge(C_j,Cell_down{5},nd1_down);
    %                     on_edge_two=on_edge(C_j,Cell_down{5},nd2_down);
    %                     on_edge_three=on_edge(C_j,Cell_down{5},nd3_down);
    %                     if on_edge_one && on_edge_two && on_edge_three
    %                         ZONE_DOWN=4;
    %                     elseif on_edge_one
    %                         ZONE_DOWN=1;
    %                     elseif on_edge_two
    %                         ZONE_DOWN=2;
    %                     elseif on_edge_three
    %                         ZONE_DOWN=3;
    %                     else
    %                         error('The zone for the stencil point is not found!');
    %                     end
    %                 end
    %                 ZONE_UP=4;
    %             else
    %                 error('The face on the boundary contains false neighbor cell info!');
    %             end
    %         else
    %             if FC{10}~=0 || FC{11}~=0
    %                 Cell_down=CELL{neigh_down(1,1)};
    %                 Cell_up=CELL{neigh_up(1,1)};
    %                 nd1_down=Cell_down{13};
    %                 nd2_down=Cell_down{14};
    %                 nd3_down=Cell_down{15};
    %                 nd1_up=Cell_up{13};
    %                 nd2_up=Cell_up{14};
    %                 nd3_up=Cell_up{15};
    %                 % down cell
    %                 in_zone_one=in_triangle(C_j_down,Cell_down{5},nd1_down,nd2_down);
    %                 in_zone_two=in_triangle(C_j_down,Cell_down{5},nd2_down,nd3_down);
    %                 in_zone_three=in_triangle(C_j_down,Cell_down{5},nd3_down,nd1_down);
    %                 if in_zone_one
    %                     ZONE_DOWN=1; % Located in zone 1
    %                 elseif in_zone_two
    %                     ZONE_DOWN=2; % Located in zone 2
    %                 elseif in_zone_three
    %                     ZONE_DOWN=3; % Located in zone 3
    %                 else
    %                     on_edge_one=on_edge(C_j_down,Cell_down{5},nd1_down);
    %                     on_edge_two=on_edge(C_j_down,Cell_down{5},nd2_down);
    %                     on_edge_three=on_edge(C_j_down,Cell_down{5},nd3_down);
    %                     if on_edge_one && on_edge_two && on_edge_three
    %                         ZONE_DOWN=4;
    %                     elseif on_edge_one
    %                         ZONE_DOWN=1;
    %                     elseif on_edge_two
    %                         ZONE_DOWN=2;
    %                     elseif on_edge_three
    %                         ZONE_DOWN=3;
    %                     else
    %                         error('The zone for the stencil point is not found!');
    %                     end
    %                 end
    %                 % up cell
    %                 in_zone_one=in_triangle(C_j_up,Cell_up{5},nd1_up,nd2_up);
    %                 in_zone_two=in_triangle(C_j_up,Cell_up{5},nd2_up,nd3_up);
    %                 in_zone_three=in_triangle(C_j_up,Cell_up{5},nd3_up,nd1_up);
    %                 if in_zone_one
    %                     ZONE_UP=1; % Located in zone 1
    %                 elseif in_zone_two
    %                     ZONE_UP=2; % Located in zone 2
    %                 elseif in_zone_three
    %                     ZONE_UP=3; % Located in zone 3
    %                 else
    %                     on_edge_one=on_edge(C_j_up,Cell_up{5},nd1_up);
    %                     on_edge_two=on_edge(C_j_up,Cell_up{5},nd2_up);
    %                     on_edge_three=on_edge(C_j_up,Cell_up{5},nd3_up);
    %                     if on_edge_one && on_edge_two && on_edge_three
    %                         ZONE_UP=4;
    %                     elseif on_edge_one
    %                         ZONE_UP=1;
    %                     elseif on_edge_two
    %                         ZONE_UP=2;
    %                     elseif on_edge_three
    %                         ZONE_UP=3;
    %                     else
    %                         error('The zone for the stencil point is not found!');
    %                     end
    %                 end
    %             end
    %         end
    %         %     %  S{10}, S{11} for V1
    %         %     DC1=S_V1{1,1};
    %         %     UC1=S_V1{2,1};
    %         %     zone_up1=zeros(1,q1);
    %         %     zone_down1=zeros(1,q1);
    %         %     for k=1:q1
    %         %         if DC1(1,k)==neigh_down(1,1) && UC1(1,k)==neigh_up(1,1)
    %         %             zone_down1(1,k)=ZONE_DOWN;
    %         %             zone_up1(1,k)=ZONE_UP;
    %         %         elseif DC1(1,k)==neigh_up(1,1) && UC1(1,k)==neigh_down(1,1)
    %         %             zone_down1(1,k)=ZONE_UP;
    %         %             zone_up1(1,k)=ZONE_DOWN;
    %         %         else
    %         %             error('Two neighbor cells for the current face reside both on one side of the face!');
    %         %         end
    %         %     end
    %         %     %  Fill the data S{10}, S{11} for V1
    %         %     S_V1{10,1}=zone_down1;
    %         %     S_V1{11,1}=zone_up1;
    %         %
    %         %     %  S{10}, S{11} for V2
    %         %     DC2=S_V2{1,1};
    %         %     UC2=S_V2{2,1};
    %         %     zone_up2=zeros(1,q2);
    %         %     zone_down2=zeros(1,q2);
    %         %     for k=1:q2
    %         %         if DC2(1,k)==neigh_down(1,1) && UC2(1,k)==neigh_up(1,1)
    %         %             zone_down2(1,k)=ZONE_DOWN;
    %         %             zone_up2(1,k)=ZONE_UP;
    %         %         elseif DC2(1,k)==neigh_up(1,1) && UC2(1,k)==neigh_down(1,1)
    %         %             zone_down2(1,k)=ZONE_UP;
    %         %             zone_up2(1,k)=ZONE_DOWN;
    %         %         else
    %         %             error('Two neighbor cells for the current face reside both on one side of the face!');
    %         %         end
    %         %     end
    %         %     %  Fill the data S{10}, S{11} for V2
    %         %     S_V2{10,1}=zone_down2;
    %         %     S_V2{11,1}=zone_up2;
    %         %% Fill in FC{16}, FC{17}
    %         FC{16}=S_V1;
    %         FC{17}=S_V2;
    %         %% FC{19}
    %         fc19=zeros(1,3);
    %         fc19(:,1)=ZONE_DOWN;
    %         fc19(:,2)=ZONE_UP;
    %         FC{19}=fc19;
    %         %% FC{20}
    %         fc20=zeros(4,3);
    %         if FC{10}==0 && FC{11}==0 % Interior face, and none of the end nodes is on boundary
    %             UP=FC{12};
    %             DOWN=FC{13};
    %             UP_Cell_number=UP(1,1);
    %             DOWN_Cell_number=DOWN(1,1);
    %             Cell_up_fixed=CELL{UP_Cell_number};
    %             C_up_fixed=Cell_up_fixed{5};
    %             Cell_down_fixed=CELL{DOWN_Cell_number};
    %             C_down_fixed=Cell_down_fixed{5};
    %             ND1=NODE{FC{8}};
    %             ND2=NODE{FC{9}};
    %             Cell_union=union(ND1{5},ND2{5});
    %             Cell_union_UP=setxor(UP_Cell_number,Cell_union);
    %             Cell_union_DOWN=setxor(DOWN_Cell_number,Cell_union);
    %             L_cell_union_up=length(Cell_union_UP);
    %             L_cell_union_down=length(Cell_union_DOWN);
    %             Pair_cell_union_up=zeros(2,L_cell_union_up*(L_cell_union_up-1)/2);
    %             Pair_cell_union_down=zeros(2,L_cell_union_down*(L_cell_union_down-1)/2);
    %             S=FC{18};
    %             C_down=S(:,1);
    %             C_up=S(:,2);
    %             %% Find the triangle that circles the stencil point in the downstream cell
    %             if single(e+norm(C_down_fixed-C_down))==single(e) % The stencil point is at the centroid
    %                 fc20(1,1)=0;
    %                 fc20(2,1)=DOWN_Cell_number;
    %                 fc20(3,1)=DOWN_Cell_number;
    %                 fc20(4,1)=DOWN_Cell_number;
    %             else
    %                 %% Find all possible paired centroids that could be used to form triangles
    %                 a_down=0;
    %                 for k=1:length(Cell_union_DOWN)
    %                     if length(Cell_union_DOWN)==2
    %                         a_down=a_down+1;
    %                         Pair_cell_union_down(:,a_down)=[Cell_union_DOWN(1);Cell_union_DOWN(2)];
    %                         break;
    %                     else
    %                         for n=2:length(Cell_union_DOWN)
    %                             a_down=a_down+1;
    %                             Pair_cell_union_down(:,a_down)=[Cell_union_DOWN(1);Cell_union_DOWN(n)];
    %                         end
    %                         Cell_union_DOWN=setxor(Cell_union_DOWN(1),Cell_union_DOWN);
    %                     end
    %                 end
    %                 %% Narrow down the previous possibility ruling out the triangles that don't
    %                 % circle the stencil point in the downstream cell
    %                 a=0;
    %                 in_tri_pair_down=0;
    %                 for k=1:a_down
    %                     Cell_pair_1=CELL{Pair_cell_union_down(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_down(2,k)};
    %                     C_down_1=Cell_pair_1{5};
    %                     C_down_2=Cell_pair_2{5};
    %                     if in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
    %                         a=a+1;
    %                         in_tri_pair_down(a)=k;
    %                     end
    %                 end
    %                 if a==0
    %                     error('No triangle is found that contains the stencil point!');
    %                 end
    %
    %                 Pair_cell_union_down_new=zeros(2,a);
    %                 for n=1:a
    %                     Pair_cell_union_down_new(:,n)=Pair_cell_union_down(:,in_tri_pair_down(n));
    %                 end
    %                 % Check
    %                 for k=1:a
    %                     Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
    %                     C_down_1=Cell_pair_1{5};
    %                     C_down_2=Cell_pair_2{5};
    %                     if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
    %                         error(':');
    %                     end
    %                 end
    %                 %% Find the pair that has the shortest distance to the stencil point
    %                 Dis_down=zeros(1,a);
    %                 for k=1:a
    %                     Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
    %                     C_down_1=Cell_pair_1{5};
    %                     C_down_2=Cell_pair_2{5};
    %                     Dis_down(1,k)=(dis(C_down,C_down_1)+dis(C_down,C_down_2))/2;
    %                 end
    %                 Dis_down_min=min(Dis_down);
    %                 for k=1:a
    %                     if single(e+Dis_down_min)==single(e+Dis_down(1,k))
    %                         break;
    %                     end
    %                 end
    %                 Cell_pair_down_found=Pair_cell_union_down_new(:,k);
    %                 % Check
    %                 Cell_pair_1=CELL{Cell_pair_down_found(1,1)};
    %                 Cell_pair_2=CELL{Cell_pair_down_found(2,1)};
    %                 C_down_1=Cell_pair_1{5};
    %                 C_down_2=Cell_pair_2{5};
    %                 if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
    %                     error(':');
    %                 end
    %                 fc20(1,1)=0;
    %                 fc20(2,1)=DOWN_Cell_number;
    %                 fc20(3,1)=Cell_pair_down_found(1,1);
    %                 fc20(4,1)=Cell_pair_down_found(2,1);
    %             end
    %             %% Find the triangle that circles the stencil point in the upstream cell
    %             if single(e+norm(C_up_fixed-C_up))==single(e) % The stencil point is at the cnetroid
    %                 fc20(1,2)=0;
    %                 fc20(2,2)=UP_Cell_number;
    %                 fc20(3,2)=UP_Cell_number;
    %                 fc20(4,2)=UP_Cell_number;
    %             else
    %                 %% Find all possible paired centroids that could be used to form triangles
    %                 a_up=0;
    %                 for k=1:length(Cell_union_UP)
    %                     if length(Cell_union_UP)==2
    %                         a_up=a_up+1;
    %                         Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(2)];
    %                         break;
    %                     else
    %                         for n=2:length(Cell_union_UP)
    %                             a_up=a_up+1;
    %                             Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(n)];
    %                         end
    %                         Cell_union_UP=setxor(Cell_union_UP(1),Cell_union_UP);
    %                     end
    %                 end
    %                 %% Narrow down the previous possibility ruling out the triangles that don't
    %                 % circle the stencil point in the upstream cell
    %                 a=0;
    %                 in_tri_pair_up=0;
    %                 for k=1:a_up
    %                     Cell_pair_1=CELL{Pair_cell_union_up(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_up(2,k)};
    %                     C_up_1=Cell_pair_1{5};
    %                     C_up_2=Cell_pair_2{5};
    %                     if in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
    %                         a=a+1;
    %                         in_tri_pair_up(a)=k;
    %                     end
    %                 end
    %                 if a==0
    %                     error('No triangle is found that contains the stencil point!');
    %                 end
    %
    %                 Pair_cell_union_up_new=zeros(2,a);
    %                 for n=1:a
    %                     Pair_cell_union_up_new(:,n)=Pair_cell_union_up(:,in_tri_pair_up(n));
    %                 end
    %                 % Check
    %                 for k=1:a
    %                     Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
    %                     C_up_1=Cell_pair_1{5};
    %                     C_up_2=Cell_pair_2{5};
    %                     if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
    %                         error(':');
    %                     end
    %                 end
    %                 %% Find the pair that has the shortest distance to the stencil point
    %                 Dis_up=zeros(1,a);
    %                 for k=1:a
    %                     Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
    %                     Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
    %                     C_up_1=Cell_pair_1{5};
    %                     C_up_2=Cell_pair_2{5};
    %                     Dis_up(1,k)=(dis(C_up,C_up_1)+dis(C_up,C_up_2))/2;
    %                 end
    %                 Dis_up_min=min(Dis_up);
    %                 for k=1:a
    %                     if single(e+Dis_up_min)==single(e+Dis_up(1,k))
    %                         break;
    %                     end
    %                 end
    %                 Cell_pair_up_found=Pair_cell_union_up_new(:,k);
    %                 % Check
    %                 Cell_pair_1=CELL{Cell_pair_up_found(1,1)};
    %                 Cell_pair_2=CELL{Cell_pair_up_found(2,1)};
    %                 C_up_1=Cell_pair_1{5};
    %                 C_up_2=Cell_pair_2{5};
    %                 if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
    %                     error(':');
    %                 end
    %                 fc20(1,2)=0;
    %                 fc20(2,2)=UP_Cell_number;
    %                 fc20(3,2)=Cell_pair_up_found(1,1);
    %                 fc20(4,2)=Cell_pair_up_found(2,1);
    %             end
    %         elseif (FC{10}~=0 && FC{11}==0) || (FC{10}==0 && FC{11}~=0) % Interior face, one of the end nodes is on boundary, the other is not
    %             % Check
    %             if FC{2}~=0
    %                 error('The current face should be interior!');
    %             end
    %             UP=FC{12};
    %             DOWN=FC{13};
    %             UP_Cell_number=UP(1,1);
    %             DOWN_Cell_number=DOWN(1,1);
    %             Cell_up_fixed=CELL{UP_Cell_number};
    %             C_up_fixed=Cell_up_fixed{5};
    %             Cell_down_fixed=CELL{DOWN_Cell_number};
    %             C_down_fixed=Cell_down_fixed{5};
    %             ND1=NODE{FC{8}};
    %             ND2=NODE{FC{9}};
    %             if FC{10}~=0 && FC{11}==0
    %                 C_nd=ND1{3};
    %                 nd=ND1{1};
    %                 Cell_union=union(ND1{5},ND2{5}); % Chould be the node on boundary
    %             elseif FC{10}==0 && FC{11}~=0
    %                 C_nd=ND2{3};
    %                 nd=ND2{1};
    %                 Cell_union=union(ND1{5},ND2{5}); % Chould be the node on boundary
    %             else
    %                 error('There must be one node on the boundary and the other is not!');
    %             end
    %             Third_cell_up=setxor(UP_Cell_number,Cell_union);
    %             Third_cell_down=setxor(DOWN_Cell_number,Cell_union);
    %             S=FC{18};
    %             C_down=S(:,1);
    %             C_up=S(:,2);
    %             %% Find the triangle that circles the stencil point in the downstream cell
    %             if single(e+norm(C_down_fixed-C_down))==single(e) % The stencil point is at the cnetroid
    %                 fc20(1,1)=0;
    %                 fc20(2,1)=DOWN_Cell_number;
    %                 fc20(3,1)=DOWN_Cell_number;
    %                 fc20(4,1)=DOWN_Cell_number;
    %             else
    %                 %% Find all possible paired centroids that could be used to form triangles
    %
    %                 %% Narrow down the previous possibility ruling out the triangles that don't
    %                 % circle the stencil point in the downstream cell
    %                 a_down=length(Third_cell_down);
    %                 a=0;
    %                 in_tri_pair_down=0;
    %                 for k=1:a_down
    %                     Cell_third=CELL{Third_cell_down(k)};
    %                     C_third_cell_down=Cell_third{5};
    %                     if in_triangle(C_down,C_down_fixed,C_nd,C_third_cell_down)
    %                         a=a+1;
    %                         in_tri_pair_down(a)=k;
    %                     end
    %                 end
    %                 if a==0
    %                     if Cell_down_fixed{10}==0 && Cell_down_fixed{11}~=0 && Cell_down_fixed{12}~=0
    %                         nd1=Cell_down_fixed{8};
    %                         nd2=Cell_down_fixed{9};
    %                         C_nd1=Cell_down_fixed{14};
    %                         C_nd2=Cell_down_fixed{15};
    %                     elseif Cell_down_fixed{10}~=0 && Cell_down_fixed{11}==0 && Cell_down_fixed{12}~=0
    %                         nd1=Cell_down_fixed{7};
    %                         nd2=Cell_down_fixed{9};
    %                         C_nd1=Cell_down_fixed{13};
    %                         C_nd2=Cell_down_fixed{15};
    %                     elseif Cell_down_fixed{10}~=0 && Cell_down_fixed{11}~=0 && Cell_down_fixed{12}==0
    %                         nd1=Cell_down_fixed{7};
    %                         nd2=Cell_down_fixed{8};
    %                         C_nd1=Cell_down_fixed{13};
    %                         C_nd2=Cell_down_fixed{14};
    %                     else
    %                         error('The down cell should be attached to boundary!');
    %                     end
    %                     if ~in_triangle(C_down,C_down_fixed,C_nd1,C_nd2)
    %                         error('The stencil point should be located in the one of the ZONE of the cell attached to boundary!');
    %                     end
    %                     fc20(1,1)=2;
    %                     fc20(2,1)=nd1;
    %                     fc20(3,1)=nd2;
    %                     fc20(4,1)=DOWN_Cell_number;
    %                 else
    %                     Third_cell_down_new=zeros(1,a);
    %                     for n=1:a
    %                         Third_cell_down_new(1,n)=Third_cell_down(in_tri_pair_down(n));
    %                     end
    %                     % Check
    %                     for k=1:a
    %                         Cell_third=CELL{Third_cell_down_new(1,k)};
    %                         C_third_cell_down=Cell_third{5};
    %                         if ~in_triangle(C_down,C_down_fixed,C_nd,C_third_cell_down)
    %                             error(':');
    %                         end
    %                     end
    %                     %% Find the pair that has the shortest distance to the stencil point
    %                     Dis_down=zeros(1,a);
    %                     for k=1:a
    %                         Cell_third=CELL{Third_cell_down_new(1,k)};
    %                         C_third_cell_down=Cell_third{5};
    %                         Dis_down(1,k)=dis(C_down,C_third_cell_down);
    %                     end
    %                     Dis_down_min=min(Dis_down);
    %                     for k=1:a
    %                         if single(e+Dis_down_min)==single(e+Dis_down(1,k))
    %                             break;
    %                         end
    %                     end
    %                     Third_cell_down_found=Third_cell_down_new(1,k);
    %                     % Check
    %                     Cell_third=CELL{Third_cell_down_found};
    %                     C_third_cell_down=Cell_third{5};
    %                     if ~in_triangle(C_down,C_down_fixed,C_nd,C_third_cell_down)
    %                         error(':');
    %                     end
    %                     fc20(1,1)=1;
    %                     fc20(2,1)=nd;
    %                     fc20(3,1)=DOWN_Cell_number;
    %                     fc20(4,1)=Third_cell_down_found;
    %                 end
    %             end
    %             %% Find the triangle that circles the stencil point in the upstream cell
    %             if single(e+norm(C_up_fixed-C_up))==single(e) % The stencil point is at the cnetroid
    %                 fc20(1,2)=0;
    %                 fc20(2,2)=UP_Cell_number;
    %                 fc20(3,2)=UP_Cell_number;
    %                 fc20(4,2)=UP_Cell_number;
    %             else
    %                 %% Find all possible paired centroids that could be used to form triangles
    %
    %                 %% Narrow down the previous possibility ruling out the triangles that don't
    %                 % circle the stencil point in the downstream cell
    %                 a_up=length(Third_cell_up);
    %                 a=0;
    %                 in_tri_pair_up=0;
    %                 for k=1:a_up
    %                     Cell_third=CELL{Third_cell_up(k)};
    %                     C_third_cell_up=Cell_third{5};
    %                     if in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
    %                         a=a+1;
    %                         in_tri_pair_up(a)=k;
    %                     end
    %                 end
    %                 if a==0
    %                     if Cell_up_fixed{10}==0 && Cell_up_fixed{11}~=0 && Cell_up_fixed{12}~=0
    %                         nd1=Cell_up_fixed{8};
    %                         nd2=Cell_up_fixed{9};
    %                         C_nd1=Cell_up_fixed{14};
    %                         C_nd2=Cell_up_fixed{15};
    %                     elseif Cell_up_fixed{10}~=0 && Cell_up_fixed{11}==0 && Cell_up_fixed{12}~=0
    %                         nd1=Cell_up_fixed{7};
    %                         nd2=Cell_up_fixed{9};
    %                         C_nd1=Cell_up_fixed{13};
    %                         C_nd2=Cell_up_fixed{15};
    %                     elseif Cell_up_fixed{10}~=0 && Cell_up_fixed{11}~=0 && Cell_up_fixed{12}==0
    %                         nd1=Cell_up_fixed{7};
    %                         nd2=Cell_up_fixed{8};
    %                         C_nd1=Cell_up_fixed{13};
    %                         C_nd2=Cell_up_fixed{14};
    %                     else
    %                         error('The up cell should be attached to boundary!');
    %                     end
    %                     if ~in_triangle(C_up,C_up_fixed,C_nd1,C_nd2)
    %                         error('The stencil point should be located in the one of the ZONE of the cell attached to boundary!');
    %                     end
    %                     fc20(1,2)=2;
    %                     fc20(2,2)=nd1;
    %                     fc20(3,2)=nd2;
    %                     fc20(4,2)=UP_Cell_number;
    %                 else
    %                     Third_cell_up_new=zeros(1,a);
    %                     for n=1:a
    %                         Third_cell_up_new(1,n)=Third_cell_up(in_tri_pair_up(n));
    %                     end
    %                     % Check
    %                     for k=1:a
    %                         Cell_third=CELL{Third_cell_up_new(1,k)};
    %                         C_third_cell_up=Cell_third{5};
    %                         if ~in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
    %                             error(':');
    %                         end
    %                     end
    %                     %% Find the pair that has the shortest distance to the stencil point
    %                     Dis_up=zeros(1,a);
    %                     for k=1:a
    %                         Cell_third=CELL{Third_cell_up_new(1,k)};
    %                         C_third_cell_up=Cell_third{5};
    %                         Dis_up(1,k)=dis(C_down,C_third_cell_up);
    %                     end
    %                     Dis_up_min=min(Dis_up);
    %                     for k=1:a
    %                         if single(e+Dis_up_min)==single(e+Dis_up(1,k))
    %                             break;
    %                         end
    %                     end
    %                     Third_cell_up_found=Third_cell_up_new(1,k);
    %                     % Check
    %                     Cell_third=CELL{Third_cell_up_found};
    %                     C_third_cell_up=Cell_third{5};
    %                     if ~in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
    %                         error(':');
    %                     end
    %                     fc20(1,2)=1;
    %                     fc20(2,2)=nd;
    %                     fc20(3,2)=UP_Cell_number;
    %                     fc20(4,2)=Third_cell_up_found;
    %                 end
    %             end
    %         elseif FC{10}~=0 && FC{11}~=0 % Boundary face
    %             % Check
    %             if FC{2}==0
    %                 error('The current face should be interior!');
    %             end
    %             fc20_bc=zeros(4,1);
    %             UP=FC{12};
    %             DOWN=FC{13};
    %             UP_Cell_number=UP(1,1);
    %             DOWN_Cell_number=DOWN(1,1);
    %             S=FC{18};
    %             if UP_Cell_number==0 && DOWN_Cell_number~=0
    %                 Cell_number=DOWN_Cell_number;
    %                 C_S=S(:,1);
    %             elseif UP_Cell_number~=0 && DOWN_Cell_number==0
    %                 Cell_number=UP_Cell_number;
    %                 C_S=S(:,2);
    %             else
    %                 error('boundary face should have only one neighbor!');
    %             end
    %             Cell_fixed=CELL{Cell_number};
    %             C_fixed=Cell_fixed{5};
    %             ND1=NODE{FC{8}};
    %             ND2=NODE{FC{9}};
    %             ND_3rd=NODE{setxor(union(FC{8},FC{9}),union(Cell_fixed{7},union(Cell_fixed{8},Cell_fixed{9})))};
    %
    %             C_nd1=ND1{3};
    %             nd1=ND1{1};
    %
    %             C_nd2=ND2{3};
    %             nd2=ND2{1};
    %
    %             Cell_union=ND_3rd{5};
    %             if single(e+norm(C_fixed-C_S))==single(e)
    %                 fc20_bc(1,1)=0;
    %                 fc20_bc(2,1)=Cell_number;
    %                 fc20_bc(3,1)=Cell_number;
    %                 fc20_bc(4,1)=Cell_number;
    %             elseif in_triangle(C_S,C_fixed,C_nd1,C_nd2)
    %                 fc20_bc(1,1)=2;
    %                 fc20_bc(2,1)=FC{8};
    %                 fc20_bc(3,1)=FC{9};
    %                 fc20_bc(4,1)=Cell_number;
    %             else
    %                 Third_cell=setxor(Cell_number,Cell_union);
    %                 b=length(Third_cell);
    %                 % Using the first end node
    %                 a1=0;
    %                 in_tri_pair_1=0;
    %                 for k=1:b
    %                     Cell_third=CELL{Third_cell(k)};
    %                     C_third_cell=Cell_third{5};
    %                     if in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
    %                         a1=a1+1;
    %                         in_tri_pair_1(a1)=k;
    %                     end
    %                 end
    %                 % Using the second end node
    %                 a2=0;
    %                 in_tri_pair_2=0;
    %                 for k=1:b
    %                     Cell_third=CELL{Third_cell(k)};
    %                     C_third_cell=Cell_third{5};
    %                     if in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
    %                         a2=a2+1;
    %                         in_tri_pair_2(a2)=k;
    %                     end
    %                 end
    %
    %                 if a1==0 && a2==0
    %                     error('There must be a triangle that can encircle the stendil point!');
    %                 elseif a1==0 && a2~=0
    %                     Third_cell_new=zeros(1,a2);
    %                     for n=1:a2
    %                         Third_cell_new(1,n)=Third_cell(in_tri_pair_2(n));
    %                     end
    %                     % Check
    %                     for k=1:a2
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
    %                             error(':');
    %                         end
    %                     end
    %                     %% Find the pair that has the shortest distance to the stencil point
    %                     Dis=zeros(1,a2);
    %                     for k=1:a2
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         Dis(1,k)=dis(C_S,C_third_cell);
    %                     end
    %                     Dis_min=min(Dis);
    %                     for k=1:a2
    %                         if single(e+Dis_min)==single(e+Dis(1,k))
    %                             break;
    %                         end
    %                     end
    %                     Third_cell_found=Third_cell_new(1,k);
    %                     % Check
    %                     Cell_third=CELL{Third_cell_found};
    %                     C_third_cell=Cell_third{5};
    %                     if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
    %                         error(':');
    %                     end
    %                     fc20_bc(1,1)=1;
    %                     fc20_bc(2,1)=nd2;
    %                     fc20_bc(3,1)=Cell_number;
    %                     fc20_bc(4,1)=Third_cell_found;
    %                 elseif a1~=0 && a2==0
    %                     Third_cell_new=zeros(1,a1);
    %                     for n=1:a1
    %                         Third_cell_new(1,n)=Third_cell(in_tri_pair_1(n));
    %                     end
    %                     % Check
    %                     for k=1:a1
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
    %                             error(':');
    %                         end
    %                     end
    %                     %% Find the pair that has the shortest distance to the stencil point
    %                     Dis=zeros(1,a1);
    %                     for k=1:a1
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         Dis(1,k)=dis(C_S,C_third_cell);
    %                     end
    %                     Dis_min=min(Dis);
    %                     for k=1:a1
    %                         if single(e+Dis_min)==single(e+Dis(1,k))
    %                             break;
    %                         end
    %                     end
    %                     Third_cell_found=Third_cell_new(1,k);
    %                     % Check
    %                     Cell_third=CELL{Third_cell_found};
    %                     C_third_cell=Cell_third{5};
    %                     if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
    %                         error(':');
    %                     end
    %                     fc20_bc(1,1)=1;
    %                     fc20_bc(2,1)=nd1;
    %                     fc20_bc(3,1)=Cell_number;
    %                     fc20_bc(4,1)=Third_cell_found;
    %                 elseif a1~=0 && a2~=0
    %                     Third_cell_new=zeros(1,a1+a2);
    %                     for n=1:a1
    %                         Third_cell_new(1,n)=Third_cell(in_tri_pair_1(n));
    %                     end
    %                     for n=1:a2
    %                         Third_cell_new(1,a1+n)=Third_cell(in_tri_pair_2(n));
    %                     end
    %                     % Check
    %                     for k=1:a1+a2
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         if k<=a1
    %                             if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
    %                                 error(':');
    %                             end
    %                         else
    %                             if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
    %                                 error(':');
    %                             end
    %                         end
    %                     end
    %                     %% Find the pair that has the shortest distance to the stencil point
    %                     Dis=zeros(1,a1+a2);
    %                     for k=1:a1+a2
    %                         Cell_third=CELL{Third_cell_new(1,k)};
    %                         C_third_cell=Cell_third{5};
    %                         Dis(1,k)=dis(C_S,C_third_cell);
    %                     end
    %                     Dis_min=min(Dis);
    %                     for k=1:a1+a2
    %                         if single(e+Dis_min)==single(e+Dis(1,k))
    %                             break;
    %                         end
    %                     end
    %                     Third_cell_found=Third_cell_new(1,k);
    %                     % Check
    %                     Cell_third=CELL{Third_cell_found};
    %                     C_third_cell=Cell_third{5};
    %                     if k<=a1
    %                         if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
    %                             error(':');
    %                         end
    %                         nd=nd1;
    %                     else
    %                         if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
    %                             error(':');
    %                         end
    %                         nd=nd2;
    %                     end
    %                     fc20_bc(1,1)=1;
    %                     fc20_bc(2,1)=nd;
    %                     fc20_bc(3,1)=Cell_number;
    %                     fc20_bc(4,1)=Third_cell_found;
    %                 else
    %                     ;
    %                 end
    %             end
    %             if UP_Cell_number==0 && DOWN_Cell_number~=0
    %                 fc20(:,1)=fc20_bc;
    %             elseif UP_Cell_number~=0 && DOWN_Cell_number==0
    %                 fc20(:,2)=fc20_bc;
    %             else
    %                 error('boundary face should have only one neighbor!');
    %             end
    %
    %         else
    %             error('Boundary identifier of end nodes on the current face is incorrect!');
    %         end
    %         FC{20}=fc20;
    %         %% FC{21}, FC{22}
    %         fc21=zeros(q1,3);
    %         fc22=zeros(q2,3);
    %         if FC{2}==0 % Interior face
    %             UP=FC{12};
    %             DOWN=FC{13};
    %             UP_Cell_number=UP(1,1);
    %             DOWN_Cell_number=DOWN(1,1);
    %             P_up=CELL{UP_Cell_number};
    %             P_down=CELL{DOWN_Cell_number};
    %             S=FC{18};
    %             Coord_down=S(:,1);
    %             Coord_up=S(:,2);
    %             v_up=Coord_up-P_up{5};
    %             v_down=Coord_down-P_down{5};
    %             %% Downwind stencil point
    %             %% V1
    %             for s=1:q1
    %                 if V1(:,s)'*v_down>=0
    %                     fc21(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %             %% V2
    %             for s=1:q2
    %                 if V2(:,s)'*v_down>=0
    %                     fc22(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %             %% upwind stencil point
    %             %% V1
    %             for s=1:q1
    %                 if V1(:,s)'*v_up>=0
    %                     fc21(s,2)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %             %% V2
    %             for s=1:q2
    %                 if V2(:,s)'*v_up>=0
    %                     fc22(s,2)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %         else % Boundary face
    %             fc21_bc=zeros(q1,1);
    %             fc22_bc=zeros(q2,1);
    %             UP=FC{12};
    %             DOWN=FC{13};
    %             UP_Cell_number=UP(1,1);
    %             DOWN_Cell_number=DOWN(1,1);
    %             S=FC{18};
    %             if UP_Cell_number==0 && DOWN_Cell_number~=0
    %                 P=CELL{DOWN_Cell_number};
    %                 Coord=S(:,1);
    %             elseif UP_Cell_number~=0 && DOWN_Cell_number==0
    %                 P=CELL{UP_Cell_number};
    %                 Coord=S(:,2);
    %             else
    %                 error('boundary face should have only one neighbor!');
    %             end
    %             v=Coord-P{5};
    %             %% V1
    %             for s=1:q1
    %                 if V1(:,s)'*v>=0
    %                     fc21_bc(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %             %% V2
    %             for s=1:q2
    %                 if V2(:,s)'*v>=0
    %                     fc22_bc(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
    %                 end
    %             end
    %             if UP_Cell_number==0 && DOWN_Cell_number~=0
    %                 fc21(:,1)=fc21_bc;
    %                 fc22(:,1)=fc22_bc;
    %             elseif UP_Cell_number~=0 && DOWN_Cell_number==0
    %                 fc21(:,2)=fc21_bc;
    %                 fc22(:,2)=fc22_bc;
    %             else
    %                 error('boundary face should have only one neighbor!');
    %             end
    %         end
    %         FC{21}=fc21;
    %         FC{22}=fc22;
    %         %% Final data filling
    %         FACE{l}=FC;
    %     end
    
    
    %% CELL type info, CL{37}
%   Type 1: the cell that none of its face is attached to inner or outer boundaries. In another word, C{2}, C{3}, and C{4} of current cell is nonzero.
%   Type 2: the cell that has one and only one face attached to the outer boundary.
%   Type 3: the cell that has one and only one face attached to the inner boundary.

    for i=1:M
        CL=CELL{i};
        if CL{2}==0 || (CL{3}==0 || CL{4}==0)
            if length(setxor([CL{2},CL{3},CL{4}],0))<2
                error(['Cell number ', num2str(i), ' has more than one face attached to boundary']);
            end
            if CL{2}==0
                FC_on_bc=FACE{CL{16}};
            elseif CL{3}==0
                FC_on_bc=FACE{CL{17}};
            elseif CL{4}==0
                FC_on_bc=FACE{CL{18}};
            else
                error('Logic error');
            end
            if FC_on_bc{24}==1 % Face on inner boundary
                CL{37}=3;
            elseif FC_on_bc{24}==2 % Face on outer boundary
                CL{37}=2;
            else
                error('logic error!');
            end
        else
            CL{37}=1;
        end
        
        CELL{i}=CL;
    end
    
    
    %% CELL info for least square calculation
    
    %% Algorithm 1 see paper Maik Stiebler, Jonas Tolke, Manfred Krafczyk. An upwind discretization scheme for the finite volume lattice. Computers & Fluids
    %% Fill CL{38~43} for cell type 1 and CL{44~50} for cell type 2 & 3
    for i=1:M
        CL=CELL{i};
        if CL{37}==1 % cell type 1
            cell1=CL{2};
            cell2=CL{3};
            cell3=CL{4};
            Cell1=CELL{cell1};
            Cell2=CELL{cell2};
            Cell3=CELL{cell3};
            dis_vec_1=Cell1{5}-CL{5};
            dis_vec_2=Cell2{5}-CL{5};
            dis_vec_3=Cell3{5}-CL{5};
            dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
            w=1./dis_sq;
            
            A=zeros(2,2);
            dis_vec=[dis_vec_1,dis_vec_2,dis_vec_3];
            for j=1:3
                A=A+w(j)*(dis_vec(:,j)*dis_vec(:,j)');
            end
            
            CL{38}=[CL{2},CL{3},CL{4}];
            CL{39}=w;
            CL{40}=dis_vec_1;
            CL{41}=dis_vec_2;
            CL{42}=dis_vec_3;
            CL{43}=A;
        else % cell type 2 and 3
            cell1=CL{2};
            cell2=CL{3};
            cell3=CL{4};
            if cell1==0
                FC_bc=FACE{CL{16}};
                Cell2=CELL{cell2};
                Cell3=CELL{cell3};
                dis_vec_1=CL{28}-CL{5};
                dis_vec_2=Cell2{5}-CL{5};
                dis_vec_3=Cell3{5}-CL{5};
            elseif cell2==0
                FC_bc=FACE{CL{17}};
                Cell1=CELL{cell1};
                Cell3=CELL{cell3};
                dis_vec_1=Cell1{5}-CL{5};
                dis_vec_2=CL{29}-CL{5};
                dis_vec_3=Cell3{5}-CL{5};
            elseif cell3==0
                FC_bc=FACE{CL{18}};
                Cell1=CELL{cell1};
                Cell2=CELL{cell2};
                dis_vec_1=Cell1{5}-CL{5};
                dis_vec_2=Cell2{5}-CL{5};
                dis_vec_3=CL{30}-CL{5};
            else
                error('Logic error!');
            end


            dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
            w=1./dis_sq;
            
            A=zeros(2,2);
            dis_vec=[dis_vec_1,dis_vec_2,dis_vec_3];
            for j=1:3
                A=A+w(j)*(dis_vec(:,j)*dis_vec(:,j)');
            end
            
            CL{44}=[CL{2},CL{3},CL{4}];
            CL{45}=[FC_bc{8},FC_bc{9}];
            CL{46}=w;
            CL{47}=dis_vec_1;
            CL{48}=dis_vec_2;
            CL{49}=dis_vec_3;
            CL{50}=A;
        end
        CELL{i}=CL;
    end
    
    %% Fill CL{51~57}
    for i=1:M
        CL=CELL{i};
        if CL{37}==2
            if CL{2}==0
                FC_on_bc=FACE{CL{16}};
            elseif CL{3}==0
                FC_on_bc=FACE{CL{17}};
            elseif CL{4}==0
                FC_on_bc=FACE{CL{18}};
            else
                error('Logic error');
            end
            if FC_on_bc{24}==2 % The face on outer boundary
                centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
                if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
                    error('logic error');
                else
                    cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
                end
                bc_cross_stencil_line=FC_on_bc{23};
                if length(setxor(bc_cross_stencil_line(2:3),0))~=1
                    error('logic error');
                else
                    bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
                end
                if bc_cross_id==1
                    centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
                elseif bc_cross_id==2
                    centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
                elseif bc_cross_id==3
                    centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
                elseif bc_cross_id==4
                    centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
                else
                    error('logic error!');
                end
                dis_vec_peri=centroid_cell_neigh_on_bc-CL{5};
                
                if CL{2}==0
                    neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
                    dis_vec_1=dis_vec_peri;
                    cell2=CL{3};
                    cell3=CL{4};
                    Cell2=CELL{cell2};
                    Cell3=CELL{cell3};
                    dis_vec_2=Cell2{5}-CL{5};
                    dis_vec_3=Cell3{5}-CL{5};
                elseif CL{3}==0
                    neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
                    dis_vec_2=dis_vec_peri;
                    cell1=CL{2};
                    cell3=CL{4};
                    Cell1=CELL{cell1};
                    Cell3=CELL{cell3};
                    dis_vec_1=Cell1{5}-CL{5};
                    dis_vec_3=Cell3{5}-CL{5};
                elseif CL{4}==0
                    neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
                    dis_vec_3=dis_vec_peri;
                    cell1=CL{2};
                    cell2=CL{3};
                    Cell1=CELL{cell1};
                    Cell2=CELL{cell2};
                    dis_vec_1=Cell1{5}-CL{5};
                    dis_vec_2=Cell2{5}-CL{5};
                else
                    error('Logic error');
                end
                
                dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
                w=1./dis_sq;
                
                A=zeros(2,2);
                dis_vec=[dis_vec_1,dis_vec_2,dis_vec_3];
                for j=1:3
                    A=A+w(j)*(dis_vec(:,j)*dis_vec(:,j)');
                end
                
                CL{51}=neigh_id;
                CL{52}=[0,0];
                CL{53}=w;
                CL{54}=dis_vec_1;
                CL{55}=dis_vec_2;
                CL{56}=dis_vec_3;
                CL{57}=A;
                
                CELL{i}=CL;
            else
                error('logic error!');
            end
        end
    end
    
    %% Fill CL{58~64}
    for i=1:M
        CL=CELL{i};
        if CL{37}==2
            CL{58}=CL{51};
            CL{59}=CL{52};
            CL{60}=CL{53};
            CL{61}=CL{54};
            CL{62}=CL{55};
            CL{63}=CL{56};
            CL{64}=CL{57};
            if CL{2}==0
                FC_on_bc=FACE{CL{16}};
            elseif CL{3}==0
                FC_on_bc=FACE{CL{17}};
            elseif CL{4}==0
                FC_on_bc=FACE{CL{18}};
            else
                error('Logic error');
            end
            if FC_on_bc{24}==2 % The face on outer boundary
                bc_cross_stencil_line=FC_on_bc{23};
                if length(setxor(bc_cross_stencil_line(2:3),0))~=1
                    error('logic error');
                else
                    bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
                end
                if bc_cross_id==1
                    CL{58}=CL{44};
                    CL{59}=CL{45};
                    CL{60}=CL{46};
                    CL{61}=CL{47};
                    CL{62}=CL{48};
                    CL{63}=CL{49};
                    CL{64}=CL{50};
                elseif bc_cross_id==2
                    ;
                elseif bc_cross_id==3
                    CL{58}=CL{44};
                    CL{59}=CL{45};
                    CL{60}=CL{46};
                    CL{61}=CL{47};
                    CL{62}=CL{48};
                    CL{63}=CL{49};
                    CL{64}=CL{50};
                elseif bc_cross_id==4
                    ;
                else
                    error('logic error!');
                end
                CELL{i}=CL;
            else
                error('logic error!');
            end
        end
    end
    
    %% Fill CL{65~71}
    for i=1:M
        CL=CELL{i};
        if CL{37}==2
            CL{65}=CL{51};
            CL{66}=CL{52};
            CL{67}=CL{53};
            CL{68}=CL{54};
            CL{69}=CL{55};
            CL{70}=CL{56};
            CL{71}=CL{57};
            if CL{2}==0
                FC_on_bc=FACE{CL{16}};
            elseif CL{3}==0
                FC_on_bc=FACE{CL{17}};
            elseif CL{4}==0
                FC_on_bc=FACE{CL{18}};
            else
                error('Logic error');
            end
            if FC_on_bc{24}==2 % The face on outer boundary
                bc_cross_stencil_line=FC_on_bc{23};
                if length(setxor(bc_cross_stencil_line(2:3),0))~=1
                    error('logic error');
                else
                    bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
                end
                if bc_cross_id==1
                    ;
                elseif bc_cross_id==2
                    CL{65}=CL{44};
                    CL{66}=CL{45};
                    CL{67}=CL{46};
                    CL{68}=CL{47};
                    CL{69}=CL{48};
                    CL{70}=CL{49};
                    CL{71}=CL{50};
                elseif bc_cross_id==3
                    ;
                elseif bc_cross_id==4
                    CL{65}=CL{44};
                    CL{66}=CL{45};
                    CL{67}=CL{46};
                    CL{68}=CL{47};
                    CL{69}=CL{48};
                    CL{70}=CL{49};
                    CL{71}=CL{50};
                else
                    error('logic error!');
                end
                CELL{i}=CL;
            else
                error('logic error!');
            end
        end
    end


      %% Algorithm 2 see paper W. Kyle Anderson and Daryl L. Bonhaus. An Implicit Upwind Algorithm for Computing Turbulent Flows on Unstructured Grids. Computers & Fluids
    %% Fill CL{38~43} for cell type 1 and CL{44~50} for cell type 2 & 3
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==1 % cell type 1
%             cell1=CL{2};
%             cell2=CL{3};
%             cell3=CL{4};
%             Cell1=CELL{cell1};
%             Cell2=CELL{cell2};
%             Cell3=CELL{cell3};
%             Coor_cent=CL{5};
%             Coor_sate=[Cell1{5},Cell2{5},Cell3{5}];
%             % r1
%             r1=0;
%             for j=1:3
%                 r1=r1+(Coor_sate(1,j)-Coor_cent(1,1))^2;
%             end
%             r1=sqrt(r1);
%             
%             % r2
%             r2=0;
%             for j=1:3
%                 r2=r2+(Coor_sate(1,j)-Coor_cent(1,1))*(Coor_sate(2,j)-Coor_cent(2,1));
%             end
%             r2=r2/r1;
%             
%             % r3
%             r3=0;
%             for j=1:3
%                 r3=r3+((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)^2;
%             end
%             r3=sqrt(r3);
%             
%             % w_x & w_y
%             w_x=zeros(1,3);
%             w_y=zeros(1,3);
%             for j=1:3
%                 w_x(1,j)=(Coor_sate(1,j)-Coor_cent(1,1))/r1^2-((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)*r2/r1/r3^2;
%                 w_y(1,j)=((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)/r3^2;
%             end
%             
%             CL{38}=[CL{2},CL{3},CL{4}];
%             CL{39}=w_x;
%             CL{40}=w_y;
%             CL{41}=0;
%             CL{42}=0;
%             CL{43}=0;
%         else % cell type 2 and 3
%             cell1=CL{2};
%             cell2=CL{3};
%             cell3=CL{4};
%             if cell1==0
%                 FC_bc=FACE{CL{16}};
%                 Cell2=CELL{cell2};
%                 Cell3=CELL{cell3};
%                 Coor_cent=CL{5};
%                 Coor_sate=[CL{28},Cell2{5},Cell3{5}];
%             elseif cell2==0
%                 FC_bc=FACE{CL{17}};
%                 Cell1=CELL{cell1};
%                 Cell3=CELL{cell3};
%                 Coor_cent=CL{5};
%                 Coor_sate=[Cell1{5},CL{29},Cell3{5}];
%             elseif cell3==0
%                 FC_bc=FACE{CL{18}};
%                 Cell1=CELL{cell1};
%                 Cell2=CELL{cell2};
%                 Coor_cent=CL{5};
%                 Coor_sate=[Cell1{5},Cell2{5},CL{30}];
%             else
%                 error('Logic error!');
%             end
%             % r1
%             r1=0;
%             for j=1:3
%                 r1=r1+(Coor_sate(1,j)-Coor_cent(1,1))^2;
%             end
%             r1=sqrt(r1);
%             
%             % r2
%             r2=0;
%             for j=1:3
%                 r2=r2+(Coor_sate(1,j)-Coor_cent(1,1))*(Coor_sate(2,j)-Coor_cent(2,1));
%             end
%             r2=r2/r1;
%             
%             % r3
%             r3=0;
%             for j=1:3
%                 r3=r3+((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)^2;
%             end
%             r3=sqrt(r3);
%             
%             % w_x & w_y
%             w_x=zeros(1,3);
%             w_y=zeros(1,3);
%             for j=1:3
%                 w_x(1,j)=(Coor_sate(1,j)-Coor_cent(1,1))/r1^2-((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)*r2/r1/r3^2;
%                 w_y(1,j)=((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)/r3^2;
%             end
%             
%             CL{44}=[CL{2},CL{3},CL{4}];
%             CL{45}=[FC_bc{8},FC_bc{9}];
%             CL{46}=w_x;
%             CL{47}=w_y;
%             CL{48}=0;
%             CL{49}=0;
%             CL{50}=0;
%         end
%         CELL{i}=CL;
%     end
%     
%     %% Fill CL{51~57}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
%                 if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
%                     error('logic error');
%                 else
%                     cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
%                 end
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
%                 elseif bc_cross_id==2
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
%                 elseif bc_cross_id==3
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
%                 elseif bc_cross_id==4
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
%                 else
%                     error('logic error!');
%                 end
%                 
%                 if CL{2}==0
%                     neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
%                     cell2=CL{3};
%                     cell3=CL{4};
%                     Cell2=CELL{cell2};
%                     Cell3=CELL{cell3};
%                     Coor_cent=CL{5};
%                     Coor_sate=[centroid_cell_neigh_on_bc,Cell2{5},Cell3{5}];
%                 elseif CL{3}==0
%                     neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
%                     cell1=CL{2};
%                     cell3=CL{4};
%                     Cell1=CELL{cell1};
%                     Cell3=CELL{cell3};
%                     Coor_cent=CL{5};
%                     Coor_sate=[Cell1{5},centroid_cell_neigh_on_bc,Cell3{5}];
%                 elseif CL{4}==0
%                     neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
%                     cell1=CL{2};
%                     cell2=CL{3};
%                     Cell1=CELL{cell1};
%                     Cell2=CELL{cell2};
%                     Coor_cent=CL{5};
%                     Coor_sate=[Cell1{5},Cell2{5},centroid_cell_neigh_on_bc];
%                 else
%                     error('Logic error');
%                 end
%                 % r1
%                 r1=0;
%                 for j=1:3
%                     r1=r1+(Coor_sate(1,j)-Coor_cent(1,1))^2;
%                 end
%                 r1=sqrt(r1);
%                 
%                 % r2
%                 r2=0;
%                 for j=1:3
%                     r2=r2+(Coor_sate(1,j)-Coor_cent(1,1))*(Coor_sate(2,j)-Coor_cent(2,1));
%                 end
%                 r2=r2/r1;
%                 
%                 % r3
%                 r3=0;
%                 for j=1:3
%                     r3=r3+((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)^2;
%                 end
%                 r3=sqrt(r3);
%                 
%                 % w_x & w_y
%                 w_x=zeros(1,3);
%                 w_y=zeros(1,3);
%                 for j=1:3
%                     w_x(1,j)=(Coor_sate(1,j)-Coor_cent(1,1))/r1^2-((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)*r2/r1/r3^2;
%                     w_y(1,j)=((Coor_sate(2,j)-Coor_cent(2,1))-(Coor_sate(1,j)-Coor_cent(1,1))*r2/r1)/r3^2;
%                 end
%                 
%                 CL{51}=neigh_id;
%                 CL{52}=[0,0];
%                 CL{53}=w_x;
%                 CL{54}=w_y;
%                 CL{55}=0;
%                 CL{56}=0;
%                 CL{57}=0;
%                 
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end
%     
%     %% Fill CL{58~64}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             CL{58}=CL{51};
%             CL{59}=CL{52};
%             CL{60}=CL{53};
%             CL{61}=CL{54};
%             CL{62}=CL{55};
%             CL{63}=CL{56};
%             CL{64}=CL{57};
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     CL{58}=CL{44};
%                     CL{59}=CL{45};
%                     CL{60}=CL{46};
%                     CL{61}=CL{47};
%                     CL{62}=CL{48};
%                     CL{63}=CL{49};
%                     CL{64}=CL{50};
%                 elseif bc_cross_id==2
%                     ;
%                 elseif bc_cross_id==3
%                     CL{58}=CL{44};
%                     CL{59}=CL{45};
%                     CL{60}=CL{46};
%                     CL{61}=CL{47};
%                     CL{62}=CL{48};
%                     CL{63}=CL{49};
%                     CL{64}=CL{50};
%                 elseif bc_cross_id==4
%                     ;
%                 else
%                     error('logic error!');
%                 end
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end
%     
%     %% Fill CL{65~71}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             CL{65}=CL{51};
%             CL{66}=CL{52};
%             CL{67}=CL{53};
%             CL{68}=CL{54};
%             CL{69}=CL{55};
%             CL{70}=CL{56};
%             CL{71}=CL{57};
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     ;
%                 elseif bc_cross_id==2
%                     CL{65}=CL{44};
%                     CL{66}=CL{45};
%                     CL{67}=CL{46};
%                     CL{68}=CL{47};
%                     CL{69}=CL{48};
%                     CL{70}=CL{49};
%                     CL{71}=CL{50};
%                 elseif bc_cross_id==3
%                     ;
%                 elseif bc_cross_id==4
%                     CL{65}=CL{44};
%                     CL{66}=CL{45};
%                     CL{67}=CL{46};
%                     CL{68}=CL{47};
%                     CL{69}=CL{48};
%                     CL{70}=CL{49};
%                     CL{71}=CL{50};
%                 else
%                     error('logic error!');
%                 end
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end
    
    %% Algorithm 3 weighted scheme based on the distance to the center centroid
    %% Fill CL{38~43} for cell type 1 and CL{44~50} for cell type 2 & 3
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==1 % cell type 1
%             cell1=CL{2};
%             cell2=CL{3};
%             cell3=CL{4};
%             Cell1=CELL{cell1};
%             Cell2=CELL{cell2};
%             Cell3=CELL{cell3};
%             dis_vec_1=Cell1{5}-CL{5};
%             dis_vec_2=Cell2{5}-CL{5};
%             dis_vec_3=Cell3{5}-CL{5};
%             dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
%             w=2/3-dis_sq/sum(dis_sq);
%             if single(sum(w))~=1
%                 error('Wrong weight');
%             end
%             
%             CL{38}=[CL{2},CL{3},CL{4}];
%             CL{39}=w;
%             CL{40}=dis_vec_1;
%             CL{41}=dis_vec_2;
%             CL{42}=dis_vec_3;
%             CL{43}=0;
%         else % cell type 2 and 3
%             cell1=CL{2};
%             cell2=CL{3};
%             cell3=CL{4};
%             if cell1==0
%                 FC_bc=FACE{CL{16}};
%                 Cell2=CELL{cell2};
%                 Cell3=CELL{cell3};
%                 dis_vec_1=CL{28}-CL{5};
%                 dis_vec_2=Cell2{5}-CL{5};
%                 dis_vec_3=Cell3{5}-CL{5};
%             elseif cell2==0
%                 FC_bc=FACE{CL{17}};
%                 Cell1=CELL{cell1};
%                 Cell3=CELL{cell3};
%                 dis_vec_1=Cell1{5}-CL{5};
%                 dis_vec_2=CL{29}-CL{5};
%                 dis_vec_3=Cell3{5}-CL{5};
%             elseif cell3==0
%                 FC_bc=FACE{CL{18}};
%                 Cell1=CELL{cell1};
%                 Cell2=CELL{cell2};
%                 dis_vec_1=Cell1{5}-CL{5};
%                 dis_vec_2=Cell2{5}-CL{5};
%                 dis_vec_3=CL{30}-CL{5};
%             else
%                 error('Logic error!');
%             end
% 
%             dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
%             w=2/3-dis_sq/sum(dis_sq);
%             if single(sum(w))~=1
%                 error('Wrong weight');
%             end
%             
%             
%             CL{44}=[CL{2},CL{3},CL{4}];
%             CL{45}=[FC_bc{8},FC_bc{9}];
%             CL{46}=w;
%             CL{47}=dis_vec_1;
%             CL{48}=dis_vec_2;
%             CL{49}=dis_vec_3;
%             CL{50}=0;
%         end
%         CELL{i}=CL;
%     end
%     
%     %% Fill CL{51~57}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 centroids_stencil_line=FC_on_bc{25}; % Change for other periodic boundary conditions
%                 if length(setxor(centroids_stencil_line(2:3),CL{1}))~=1
%                     error('logic error');
%                 else
%                     cell_neigh_on_bc=CELL{setxor(centroids_stencil_line(2:3),CL{1})};
%                 end
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[0;(Y2-Y1)];
%                 elseif bc_cross_id==2
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}+[(X2-X1);0];
%                 elseif bc_cross_id==3
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[0;(Y2-Y1)];
%                 elseif bc_cross_id==4
%                     centroid_cell_neigh_on_bc=cell_neigh_on_bc{5}-[(X2-X1);0];
%                 else
%                     error('logic error!');
%                 end
%                 dis_vec_peri=centroid_cell_neigh_on_bc-CL{5};
%                 
%                 if CL{2}==0
%                     neigh_id=[cell_neigh_on_bc{1},CL{3},CL{4}];
%                     dis_vec_1=dis_vec_peri;
%                     cell2=CL{3};
%                     cell3=CL{4};
%                     Cell2=CELL{cell2};
%                     Cell3=CELL{cell3};
%                     dis_vec_2=Cell2{5}-CL{5};
%                     dis_vec_3=Cell3{5}-CL{5};
%                 elseif CL{3}==0
%                     neigh_id=[CL{2},cell_neigh_on_bc{1},CL{4}];
%                     dis_vec_2=dis_vec_peri;
%                     cell1=CL{2};
%                     cell3=CL{4};
%                     Cell1=CELL{cell1};
%                     Cell3=CELL{cell3};
%                     dis_vec_1=Cell1{5}-CL{5};
%                     dis_vec_3=Cell3{5}-CL{5};
%                 elseif CL{4}==0
%                     neigh_id=[CL{2},CL{3},cell_neigh_on_bc{1}];
%                     dis_vec_3=dis_vec_peri;
%                     cell1=CL{2};
%                     cell2=CL{3};
%                     Cell1=CELL{cell1};
%                     Cell2=CELL{cell2};
%                     dis_vec_1=Cell1{5}-CL{5};
%                     dis_vec_2=Cell2{5}-CL{5};
%                 else
%                     error('Logic error');
%                 end
%                 
%                 dis_sq=[dis_vec_1'*dis_vec_1,dis_vec_2'*dis_vec_2,dis_vec_3'*dis_vec_3];
%                 w=2/3-dis_sq/sum(dis_sq);
%                 if single(sum(w))~=1
%                     error('Wrong weight');
%                 end
%                 
%                 CL{51}=neigh_id;
%                 CL{52}=[0,0];
%                 CL{53}=w;
%                 CL{54}=dis_vec_1;
%                 CL{55}=dis_vec_2;
%                 CL{56}=dis_vec_3;
%                 CL{57}=0;
%                 
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end
%     
%     %% Fill CL{58~64}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             CL{58}=CL{51};
%             CL{59}=CL{52};
%             CL{60}=CL{53};
%             CL{61}=CL{54};
%             CL{62}=CL{55};
%             CL{63}=CL{56};
%             CL{64}=CL{57};
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     CL{58}=CL{44};
%                     CL{59}=CL{45};
%                     CL{60}=CL{46};
%                     CL{61}=CL{47};
%                     CL{62}=CL{48};
%                     CL{63}=CL{49};
%                     CL{64}=CL{50};
%                 elseif bc_cross_id==2
%                     ;
%                 elseif bc_cross_id==3
%                     CL{58}=CL{44};
%                     CL{59}=CL{45};
%                     CL{60}=CL{46};
%                     CL{61}=CL{47};
%                     CL{62}=CL{48};
%                     CL{63}=CL{49};
%                     CL{64}=CL{50};
%                 elseif bc_cross_id==4
%                     ;
%                 else
%                     error('logic error!');
%                 end
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end
%     
%     %% Fill CL{65~71}
%     for i=1:M
%         CL=CELL{i};
%         if CL{37}==2
%             CL{65}=CL{51};
%             CL{66}=CL{52};
%             CL{67}=CL{53};
%             CL{68}=CL{54};
%             CL{69}=CL{55};
%             CL{70}=CL{56};
%             CL{71}=CL{57};
%             if CL{2}==0
%                 FC_on_bc=FACE{CL{16}};
%             elseif CL{3}==0
%                 FC_on_bc=FACE{CL{17}};
%             elseif CL{4}==0
%                 FC_on_bc=FACE{CL{18}};
%             else
%                 error('Logic error');
%             end
%             if FC_on_bc{24}==2 % The face on outer boundary
%                 bc_cross_stencil_line=FC_on_bc{23};
%                 if length(setxor(bc_cross_stencil_line(2:3),0))~=1
%                     error('logic error');
%                 else
%                     bc_cross_id=setxor(bc_cross_stencil_line(2:3),0);
%                 end
%                 if bc_cross_id==1
%                     ;
%                 elseif bc_cross_id==2
%                     CL{65}=CL{44};
%                     CL{66}=CL{45};
%                     CL{67}=CL{46};
%                     CL{68}=CL{47};
%                     CL{69}=CL{48};
%                     CL{70}=CL{49};
%                     CL{71}=CL{50};
%                 elseif bc_cross_id==3
%                     ;
%                 elseif bc_cross_id==4
%                     CL{65}=CL{44};
%                     CL{66}=CL{45};
%                     CL{67}=CL{46};
%                     CL{68}=CL{47};
%                     CL{69}=CL{48};
%                     CL{70}=CL{49};
%                     CL{71}=CL{50};
%                 else
%                     error('logic error!');
%                 end
%                 CELL{i}=CL;
%             else
%                 error('logic error!');
%             end
%         end
%     end