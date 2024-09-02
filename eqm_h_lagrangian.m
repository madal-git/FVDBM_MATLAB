function f_eq_lagrangian = eqm_h_lagrangian(CELL,M,f,f_node,f_eq,V,V1,V2,q,w,Rho_r,fmp,fpdc,fd)
% f_eq_lagrandian = eqm_h_lagrangian (CL,f,V,q,w,Rho_r,fd) calculates the
% hydrodynamic equilibrium pdf at the cell barycenters of all cells in Lagrangian fashion
% according to given lattice
% CELL is the entire CELL data structure
% M is the total number of cells
% f is the entire hydrodynamic pdf data at centroids
% f_node is the entire hydrodynamic pdf data at nodes
% f_eq is the entire hydrodynamic quilibrium pdf data at nodes
% V is the currently applied velocity decomposition of the lattice
% V1 is the fisrt possible velocity decomposition of the lattice
% V2 is the second possible velocity decomposition of the lattice
% q is the total number of velocity components in the lattice
% w is the weighting factor for calculating equilibrium pdf accoring to 
% the lattice structure
% Rho_f is the reference density
% fmp is the flag for different mapping method
% fpdc is the flag for periodic boundary conditions. fpdc=0---No periodic
% boundaries; fpdc=1---Only left & right boundaries are periodic;
% fpdc=2---Only top & bottom boundaries are periodic; fpdc=3---All
% boundaries are periodic
% fd is the flag for which density is used. fd=0----local density; fd=1----reference density

if sum(sum(V-V1))==0
    cell_data_pointer=57;
elseif sum(sum(V-V2))==0
    cell_data_pointer=58;
else
    error('Wrong lattice!');
end

f_lagrangian=zeros(q,M);
f_eq_lagrangian=f_lagrangian;

%% Calculate Lagrangian pdf for each cell
% A=zeros(2,M);
fmp=1;
if fmp==1
    f_lagrangian=f;
else
    if fmp==6
        gradient=0; % Change for fmp=6
    else
        gradient=0;
    end
    for k=1:M
        CL=CELL{k};
        coord_lagrangian=CL{cell_data_pointer};
        for j=1:q
%             if j==1
%                 [f_lagrandian(j,k),A(:,k)]=in_cell_mapping_temp('dummy',CELL,'dummy',f(j,:),f_node(j,:),k,coord_lagrangian(:,j),'dummy','dummy',gradient,fmp,V,fpdc);
%             else
            f_lagrangian(j,k)=in_cell_mapping('dummy',CELL,'dummy',f(j,:),f_node(j,:),k,coord_lagrangian(:,j),'dummy','dummy',gradient,fmp,V,fpdc);
%             end
        end
    end
end


if fmp==1
    f_eq_lagrangian=f_eq;
else
    %% Calculate macros for each cell
    [RHO,U,T]=macro_h(f_lagrangian,V,Rho_r,fd);
    %% Calculate the Lagrangian f_eq for implicit f_eq at t_n+1
    for k=1:M
        f_eq_lagrangian(:,k)=eqm_h(V,U(:,k),RHO(1,k),T(1,k),q,w,Rho_r,fd);
    end
end
