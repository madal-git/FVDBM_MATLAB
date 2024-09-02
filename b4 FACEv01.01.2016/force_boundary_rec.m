function [F_p_total,F_m_n_total,F_m_s_total]=force_boundary_rec(boundary,X1,X2,Y1,Y2,NODE,CELL,FACE,Rho,U,Rho_nd,U_nd,V,Tau,controller)
% function [F_p_total,F_m_n_total,F_m_s_total]=force_boundary_rec(boundary,NODE,CELL,FACE,Rho,U,Rho_nd,U_nd,V,Tau,controller)
% returns the total force on any given boundary of a rectangular domain
% F_p_total is the force due to pressure on that boundary, a column vector [x;y]
% F_m_n_total is the force due to the fluid velocity gradient normal to the
% boundary surface(normal stress), a column vector [x;y]
% F_m_s_total is the force due to the fluid velocity gradient tangantial to
% the boundary surface(shear stress), a column vector [x;y]
% boundary is the name of given boundary, it has to be Top, Right, Bottom,
% Left or Immersed
% X1, X2, Y1 and Y2 are the boundary limits of the rectangular domain.
% NODE is the NODE data structure
% CELL is the CELL triangle data structure
% FACE is the FACE data structure
% Rho is the density data of all triangles
% U is the velocity data of all triangles
% Rho_nd is the density data of all nodes
% U_nd is the velocity data of all nodes
% V is the lattice applied.
% Tau is the re;axation time;
% controller is the flag for differet scheme of calculating force.
% 0---Mathematical approach; 1---Physical approach; 2---Hybrid

%% Check the legitimacy of given boundary name and prepare the edge and tri data for function force_edge
e=100*(X2+Y2)/2;
M=length(CELL);
if strcmp(boundary,'Top')==1
    edge_found_counter=0;
    for r=1:M
        P=CELL{r};
        boundary_edge_found=0;
        for i=1:3
            Mid=P{28+i-1};
            if single(e+Mid(2,1))==single(e+Y2)
                boundary_edge_found=1;
                break;
            end
        end
        if boundary_edge_found==1
            edge_found_counter=edge_found_counter+1;
            tri(edge_found_counter)=r;
            edge(edge_found_counter)=i;
        end
    end
elseif strcmp(boundary,'Right')==1
    edge_found_counter=0;
    for r=1:M
        P=CELL{r};
        boundary_edge_found=0;
        for i=1:3
            Mid=P{28+i-1};
            if single(e+Mid(1,1))==single(e+X2)
                boundary_edge_found=1;
                break;
            end
        end
        if boundary_edge_found==1
            edge_found_counter=edge_found_counter+1;
            tri(edge_found_counter)=r;
            edge(edge_found_counter)=i;
        end
    end
elseif strcmp(boundary,'Bottom')==1
    edge_found_counter=0;
    for r=1:M
        P=CELL{r};
        boundary_edge_found=0;
        for i=1:3
            Mid=P{28+i-1};
            if single(e+Mid(2,1))==single(e+Y1)
                boundary_edge_found=1;
                break;
            end
        end
        if boundary_edge_found==1
            edge_found_counter=edge_found_counter+1;
            tri(edge_found_counter)=r;
            edge(edge_found_counter)=i;
        end
    end
elseif strcmp(boundary,'Left')==1
    edge_found_counter=0;
    for r=1:M
        P=CELL{r};
        boundary_edge_found=0;
        for i=1:3
            Mid=P{28+i-1};
            if single(e+Mid(1,1))==single(e+X1)
                boundary_edge_found=1;
                break;
            end
        end
        if boundary_edge_found==1
            edge_found_counter=edge_found_counter+1;
            tri(edge_found_counter)=r;
            edge(edge_found_counter)=i;
        end
    end
elseif strcmp(boundary,'Immersed')==1
    edge_found_counter=0;
    for r=1:M
        P=CELL{r};
        boundary_edge_found=0;
        for i=1:3
            if single(P{19+(i-1)})<0
                boundary_edge_found=1;
                break;
            end
        end
        if boundary_edge_found==1
            edge_found_counter=edge_found_counter+1;
            tri(edge_found_counter)=r;
            edge(edge_found_counter)=i;
        end
    end
else
    error('Please check the spelling of the boundary name!');
end

%% Calculate the total force on the given boundary
% Check tri and edge
if length(tri)~=edge_found_counter || length(edge)~=edge_found_counter
    error('The numbder of found edges is not correct!');
end
% calculate force
F_p_total=[0;0];
F_m_n_total=[0;0];
F_m_s_total=[0;0];
for s=1:edge_found_counter
    [F_p,F_m_n,F_m_s]=force_edge(edge(s),tri(s),NODE,CELL,FACE,Rho,U,Rho_nd,U_nd,V,Tau,controller);
    F_p_total=F_p_total+F_p;
    F_m_n_total=F_m_n_total+F_m_n;
    F_m_s_total=F_m_s_total+F_m_s;
end
