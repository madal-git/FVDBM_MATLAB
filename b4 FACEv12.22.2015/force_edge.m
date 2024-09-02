function [F_p,F_m_n,F_m_s]=force_edge(edge,tri,NODE,CELL,FACE,Rho,U,Rho_nd,U_nd,V,Tau,controller)
% function [F_p,F_m_n,F_m_s]=force_edge(nd1,nd2,r,NODE,CELL,UWD,controller)
% returns the force evaluations on the given single edge
% F_p is the force due to pressure on that edge, a column vector [x;y]
% F_m_n is the force due to the fluid velocity gradient normal to the edge
% (normal stress), a column vector [x;y]
% F_m_s is the force due to the fluid velocity gradient tangantial to the
% edge (shear stress), a column vector [x;y]
% edge is the edge-th edge of triangle tri (1,2 or 3)
% tri is the triangle number that contains nd1 and nd2
% NODE is the NODE data structure
% CELL is the cell data structure
% FACE is the FACE data structure
% Rho is the density data of all triangles
% U is the velocity data of all triangles
% Rho_nd is the density data of all nodes
% U_nd is the velocity data of all nodes
% V is the lattice applied.
% Tau is the re;axation time;
% controller is the flag for differet scheme of calculating force.
% 0---Mathematical approach; 1---Physical approach

%% Check whether the edge is on boundary
P=CELL{tri};
if P{19+edge-1}==0
    error('The current edge is not on boundary, cannot evaluate force!');
end
%% Force calculation
% determine the two nodes at the ends of the edge
if edge==1
    ND1=NODE{P{7}};
    ND2=NODE{P{8}};
elseif edge==2
    ND1=NODE{P{8}};
    ND2=NODE{P{9}};
elseif edge==3
    ND1=NODE{P{9}};
    ND2=NODE{P{7}};
else
    error('The given edge number is incorrect!');
end
% the macro varibles on the edge
Rho_edge=(Rho_nd(ND1{1})+Rho_nd(ND2{1}))/2;
U_edge=(U_nd(:,ND1{1})+U_nd(:,ND2{1}))/2;
%  n_U_edge=U_edge/sqrt(U_edge(1)^2+U_edge(2)^2);
% the macro varibles that on the perpendicular line that cross the center
% point of the edge (The same assumption as that for flux calculation, which assumes
% the value at the projected point is equal to the centroid)
Rho_p=Rho(tri);
U_p=U(:,tri);
%  n_U_p=U_p/sqrt(U_p(1)^2+U_p(2)^2);
% geometry of the edge
l_edge=P{31+edge-1};
n_edge_n=-P{34+edge-1}';
n_edge_t=[n_edge_n(2);-n_edge_n(1)];
dis_edge2c=dis(P{28+edge-1},P{5});
n_edge2c=(P{5}-P{28+edge-1})/dis_edge2c;
UD=FACE{P{16+edge-1}};
US=UD{16};
if mean(US{17})~=2
    error('The current edge is not on boundary');
end
dis_edge2p=mean(US{13})/2;
% Calculation executed
q=length(V(1,:));
if q==7
    F_p=-Rho_edge/4*n_edge_n*l_edge;
    Mew=Tau/4;
elseif q==9
    F_p=-Rho_edge/3*n_edge_n*l_edge;
    Mew=Tau/3;
elseif q==13
    F_p=-Rho_edge/2*n_edge_n*l_edge;
    Mew=Tau/2;
else
    error('Other lattice is not available!');
end
if controller==0
    % normal
    dUn_over_dn_edge2c=((U_p(1)*n_edge_n(1)+U_p(2)*n_edge_n(2))-(U_edge(1)*n_edge_n(1)+U_edge(2)*n_edge_n(2)))/dis_edge2c;
    dUn_over_dn=dUn_over_dn_edge2c*(n_edge_n(1)*n_edge2c(1)+n_edge_n(2)*n_edge2c(2));
    F_m_n=dUn_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_n;
    % tangential
    dUt_over_dn_edge2c=((U_p(1)*n_edge_t(1)+U_p(2)*n_edge_t(2))-(U_edge(1)*n_edge_t(1)+U_edge(2)*n_edge_t(2)))/dis_edge2c;
    dUt_over_dn=dUt_over_dn_edge2c*(n_edge_n(1)*n_edge2c(1)+n_edge_n(2)*n_edge2c(2));
    F_m_s=dUt_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_t;
elseif controller==1
    % normal
    dUn_over_dn=((U_p(1)*n_edge_n(1)+U_p(2)*n_edge_n(2))-(U_edge(1)*n_edge_n(1)+U_edge(2)*n_edge_n(2)))/dis_edge2p;
    F_m_n=dUn_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_n;
    % tangential
    dUt_over_dn=((U_p(1)*n_edge_t(1)+U_p(2)*n_edge_t(2))-(U_edge(1)*n_edge_t(1)+U_edge(2)*n_edge_t(2)))/dis_edge2p;
    F_m_s=dUt_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_t;
elseif controller==2 % hybrid of 1 and 2
    % normal
    dUn_over_dn_edge2c=((U_p(1)*n_edge_n(1)+U_p(2)*n_edge_n(2))-(U_edge(1)*n_edge_n(1)+U_edge(2)*n_edge_n(2)))/dis_edge2c;
    dUn_over_dn=dUn_over_dn_edge2c*(n_edge_n(1)*n_edge2c(1)+n_edge_n(2)*n_edge2c(2));
    F_m_n=dUn_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_n;
    % tangential
    dUt_over_dn=((U_p(1)*n_edge_t(1)+U_p(2)*n_edge_t(2))-(U_edge(1)*n_edge_t(1)+U_edge(2)*n_edge_t(2)))/dis_edge2p;
    F_m_s=dUt_over_dn*(Rho_edge+Rho_p)/2*Mew*l_edge*n_edge_t;
else
    error('Incorrect controller for force evauation!');
end
