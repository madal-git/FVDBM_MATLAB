function fl=flux_pdf_t(CELL,M,NODE,N,FACE,O,g_nd,g_old,geq,RHO,Rho_r,V,dt,wl,~,wt,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,fppi,ftvd,fpdc,fd,fm,fmp)
% fl=flux_pdf_t(CELL,M,NODE,N,FACE,O,g_nd,g_old,geq,RHO,Rho_r,V,dt,wl,wlb,wt,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,fppi,ftvd,fpdc,fd,fm,fmp)
% calculates the total flux of pdf occured to each trianglar cell.
% CELL is the cell data structure
% M is the length of CELL
% NODE is the node data structure
% N is the length of NODE
% FACE is the face data structure
% O is the length of FACE
% g_nd is the up-to-date entire nodal thermal pdf matrix
% g_old is the up-to-date entire cell centroid pdf matrix
% geq is the up-to-date entire cell equilibrium pdf matrix
% RHO is the up-to-date entire cell equilibrium pdf vector
% Rho_r is the reference density
% V is the lattice actually applied
% wl is the weighting factor of average nodal flux, then (1-wl) is the
% weight of upwind flux
% wlb is the weighting factor of average nodal flux on boundary, then (1-wlb) is the
% weight of upwind flux on boundary
% V1 is the first lattice structure in stencil.m
% V2 is the second lattice structure in stencil.m

% interior is one-digit numeric flag for different scheme of flux through internal edges
% periodic is 1-digit numeric flag for different scheme of flux through Periodic edges
% inlet is 1-digit numeric flag for different scheme of flux through Velocity Inlet edges
% outlet is 1-digit numeric flag for different scheme of flux through Velocity Outlet edges
% s_wall is 1-digit numeric flag for different scheme of flux through Stationary Wall edges
% m_wall is 1-digit numeric flag for different scheme of flux through Moving Wall edges
% well_dp is 1-digit numeric flag for different scheme of flux through Zero Gradient edges
% 0---1st-order upwind (FOU)
% 1---Godunov-PL, Lax-Wendroff
% 2---Godunov-PL, Beam-Warming
% 3---Godunov-PL, Fromm
% 4---Non-Godunov general form
% 5---Non-Godunov,2nd-order upwind (SOU)
% 6---Non-Godunov, TVD
% 7---Non-Godunov,QUICK
% 8---Non-Godunov,QUICKEST
% 9---Godunov-PP
% 10--Godunov, TVD, PC + PL-LW
% 11--Godunov, TVD, PC + PP
% 12--Godunov, TVD, PL-LW + PP

% FPPI is the coefficient for the interpolation of PDF on the face between
% upwind anf further upwind cells. FPPI is positve integer ONLY

% FTVD is the flag for TVD scheme

% FPDC is the flag for periodic boundary conditions. FPDC=0---No periodic
% boundaries; FPDC=1---Only left & right boundaries are periodic;
% FPDC=2---Only top & bottom boundaries are periodic; FPDC=3---All
% boundaries are periodic

% fd is the flag for which density is used in eqm_h. fd=0----local density; fd=1----reference density

% fm is the flag for mesh type. fm=0---IRT mesh; fm=1---general mesh
% fmp is the flag for mapping
% FMP=1---Global 1st-order mapping
% FMP=2---Global 2nd-order mapping
% FMP=3---1st-order mapping for the stencil point with two enclosing boundary nodes, 2nd-order mapping for other stencil points
% FMP=4---1st-order mapping for the stecil point with one or two enclosing boundary nodes, 2nd-order mapping for other stencil points
% FMP=5---1st-order mapping for the stencil point with all three centroid enclosing points, 2nd-order mapping for other stencil points

% fl is the matrix of total flux of pdf through the boundaries of each
% trianglar CV for all cells.
% fl=wl*nfl+(1-wl)*ufl

q=length(g_nd(:,1));
if q~=length(V(1,:))
    error('The decomposition of pdf is not compatible with the lattice structure');
end

flux_cell_marker=zeros(3,M);
fl=zeros(q,M);
fl_face=zeros(q,O);

%% Determing the lattice applied
if norm(V-V1,1)==0
    lattice=1;
elseif norm(V-V2,1)==0
    lattice=2;
else
    error('The lattice used is not the lattice used for determing upwind cells');
end

for r=1:O
    %% Basic info
    FC=FACE{r};
    face_bc_flag=FC{2};
    face_type=FC{24};
    face_normal_lattice=FC{4}*V;
    %% Choose lattice-dependent info
    if lattice==1
        ls_no_pdc=FC{46};
        ls_all_pdc=FC{48};
        ls_lr_pdc=FC{50};
        ls_tb_pdc=FC{52};
        latt_dir_marker=FC{54};
    elseif lattice==2
        ls_no_pdc=FC{47};
        ls_all_pdc=FC{49};
        ls_lr_pdc=FC{51};
        ls_tb_pdc=FC{53};
        latt_dir_marker=FC{55};
    else
        error('The lattice used is not the lattice used for determing upwind cells');
    end
    
    %% Select stencil and lattice stencil info based on the periodic conditions
    if fpdc==0 % No periodic boundaries
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=FC{21};
        
        s_cell_234=FC{16};
        s_stcl_coord_234=FC{17};
        s_b_intcp_234=FC{18};
        s_map_id_234=FC{20};
        s_map_coord_234=FC{21};
        
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_no_pdc;
    elseif fpdc==1 % Left and right boundaries are periodic
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=FC{21};
        
        s_cell_234=FC{31};
        s_stcl_coord_234=FC{32};
        s_b_intcp_234=FC{34};
        s_map_id_234=FC{36};
        s_map_coord_234=FC{37};
        
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_lr_pdc;
    elseif fpdc==2 % Top and bottom boundaries are periodic
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=FC{21};
        
        s_cell_234=FC{38};
        s_stcl_coord_234=FC{39};
        s_b_intcp_234=FC{41};
        s_map_id_234=FC{43};
        s_map_coord_234=FC{44};
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_tb_pdc;
    elseif fpdc==3 % All outer boundaries are periodic
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=FC{21};
        
        s_cell_234=FC{25};
        s_stcl_coord_234=FC{26};
        s_b_intcp_234=zeros(1,length(FC{16}));
        s_map_id_234=FC{29};
        s_map_coord_234=FC{30};
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_all_pdc;
    else
        error('Wrong flag for periodic boundary conditions!');
    end
    
    if face_bc_flag==0 % Interior face
        %% Calculate the PDFs at each stencil points
        if interior<=1 % 1st-order upwind and LW
            if face_type==1
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
            elseif face_type==3 || face_type==4
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
            else
                error('Wrong flag for interior faces!');
            end
            fs_dd=fs_d; % Dummy
            fs_uu=fs_u; % Dummy
        else % Other flux calculation schemes that requires all four stencil points
            if face_type==1
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                if s_cell_1(:,1)==0
                    fs_intcpt=(1-s_b_intcp_1(3,1))*g_nd(:,s_b_intcp_1(1,1))+s_b_intcp_1(3,1)*g_nd(:,s_b_intcp_1(2,1));
                    fs_dd=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),fmp,V);
                end
                if s_cell_1(:,4)==0
                    fs_intcpt=(1-s_b_intcp_1(3,4))*g_nd(:,s_b_intcp_1(1,4))+s_b_intcp_1(3,4)*g_nd(:,s_b_intcp_1(2,4));
                    fs_uu=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                else
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),fmp,V);
                end
            elseif face_type==3 || face_type==4
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                if s_cell_234(:,1)==0
                    fs_intcpt=(1-s_b_intcp_234(3,1))*g_nd(:,s_b_intcp_234(1,1))+s_b_intcp_234(3,1)*g_nd(:,s_b_intcp_234(2,1));
                    fs_dd=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                end
                if s_cell_234(:,4)==0
                    fs_intcpt=(1-s_b_intcp_234(3,4))*g_nd(:,s_b_intcp_234(1,4))+s_b_intcp_234(3,4)*g_nd(:,s_b_intcp_234(2,4));
                    fs_uu=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                else
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                end
            else
                error('Wrong flag for interior faces!');
            end
        end
%         f_stencil=[fs_dd,fs_d,fs_u,fs_uu];
        %% calculate the PDFs at the face center
        if face_type==1 || face_type==4
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_14,latt_dir_marker,dt,face_normal_lattice,interior,fppi,ftvd);
        elseif face_type==3
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,interior,fppi,ftvd);
        else
            error('Logic error!');
        end
    elseif face_bc_flag>0 % Outer boundary face
        if face_bc_flag==1 % Periodic
            %% Calculate the PDFs at each stencil points
            if periodic<=1 % 1st-order upwind and LW
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,periodic,fppi,ftvd);
        elseif face_bc_flag==20 || face_bc_flag==21 % Velocity inlet or pressure inlet
            %% Calculate the PDFs at each stencil points
            if inlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,inlet,fppi,ftvd);
        elseif face_bc_flag==30 || face_bc_flag==31 % Velocity outlet or pressure outlet
            %% Calculate the PDFs at each stencil points
            if outlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,outlet,fppi,ftvd);
        elseif face_bc_flag==4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
%             f_stencil=[fs_dd,fs_d,fs_u,fs_uu];
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,s_wall,fppi,ftvd);
        elseif face_bc_flag==5
            %% Calculate the PDFs at each stencil points
            if m_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*g_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*g_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*g_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*g_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,m_wall,fppi,ftvd);
        elseif face_bc_flag==6
            %% Calculate the PDFs at each stencil points
            if well_dp<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                    fs_dd=fs_d; % Zero gradient
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                    fs_uu=fs_u; % Zero gradient
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,well_dp,fppi,ftvd);
        else
            error('Wrong flag for face on outer boundaries!');
        end
    elseif face_bc_flag<0 % Inner boundary face
        if face_bc_flag==-4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*g_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*g_nd(:,s_b_intcp_1(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*g_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*g_nd(:,s_b_intcp_1(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*g_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*g_nd(:,s_b_intcp_1(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,g_old,g_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*g_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*g_nd(:,s_b_intcp_1(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_14,latt_dir_marker,dt,face_normal_lattice,s_wall,fppi,ftvd);
        else
            error('Wrong flag for face on outer boundaries!');
        end
    else
        error('Logic Error!');
    end
    fl_face(:,r)=(f_face.*(face_normal_lattice)')*FC{3};
    %fl_face(:,r)=diag(f_face*(face_normal_lattice))*FC{3};
end

%%%% Deliver flux on each to each cell
for r=1:O
    FC=FACE{r};
    cell_1=FC{12};
    cell_2=FC{13};
    %%%%% First neighbor cell
    if cell_1(1,1)~=0
        C1=CELL{cell_1(1,1)};
        n1=C1{33+cell_1(1,2)};
        bool1=FC{4}*n1';
        if single(100+bool1)==single(100-1)
            fl(:,cell_1(1,1))=fl(:,cell_1(1,1))-fl_face(:,r);
        elseif single(100+bool1)==single(100+1)
            fl(:,cell_1(1,1))=fl(:,cell_1(1,1))+fl_face(:,r);
        else
            error('The face unit normal is incorrect!');
        end
        %%%%% Mark marker
        flux_cell_marker(cell_1(1,2),cell_1(1,1))=flux_cell_marker(cell_1(1,2),cell_1(1,1))+1;
    end
    if cell_2(1,1)~=0
        C2=CELL{cell_2(1,1)};
        n2=C2{33+cell_2(1,2)};
        bool2=FC{4}*n2';
        %%%%% Second neighbor cell
        if single(100+bool2)==single(100-1)
            fl(:,cell_2(1,1))=fl(:,cell_2(1,1))-fl_face(:,r);
        elseif single(100+bool2)==single(100+1)
            fl(:,cell_2(1,1))=fl(:,cell_2(1,1))+fl_face(:,r);
        else
            error('The face unit normal is incorrect!');
        end
        flux_cell_marker(cell_2(1,2),cell_2(1,1))=flux_cell_marker(cell_2(1,2),cell_2(1,1))+1;
    end
end
%%%% Check
if mean(mean(flux_cell_marker))~=1
    error('flux is missing on some faces!');
end