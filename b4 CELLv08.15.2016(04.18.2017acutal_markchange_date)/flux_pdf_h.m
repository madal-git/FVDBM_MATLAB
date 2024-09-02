function fl=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,feq,RHO,Rho_r,V,dt,wl,wlb,wh,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,iof,fupd,ftvd,fpdc,fd,fm,fmp,fe)
% fl=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,feq,RHO,Rho_r,V,dt,wl,wlb,wh,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,fupd,ftvd,fpdc,fd,fm,fmp,fe)
% calculates the total flux of pdf occured to each trianglar cell.
% CELL is the cell data structure
% M is the length of CELL
% NODE is the node data structure
% N is the length of NODE
% FACE is the face data structure
% O is the length of FACE
% f_nd is the up-to-date entire nodal pdf matrix
% f_old is the up-to-date entire cell centroid pdf matrix
% feq is the up-to-date entire cell equilibrium pdf matrix
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

% FUPD is the flag for wether use upwind scheme to calculate the flux

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
% FMP=6---least-square

% fe is the flag for extrapolation of pdf along the stencil for the ghost
% stencil points on boundary
% fe=0---linear extrapolation
% fe=1---cubic extrapolation

% fl is the matrix of total flux of pdf through the boundaries of each
% trianglar CV for all cells.
% fl=wl*nfl+(1-wl)*ufl

q=length(f_nd(:,1));
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

%% Calculate the gradient of PDF at cell centroids
if fmp==6
    gf=gradient_tri(q,M,CELL,f_old,f_nd,fpdc);
else
    gf=0;
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
    s_dis=FC{22};
    
    if face_bc_flag==0 % Interior face
        %% Calculate the PDFs at each stencil points
        if interior<=1 % 1st-order upwind and LW
            if face_type==1
                %% Create the pdf values at the stencil points
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),gf,fmp,V);
                %% Create the pdf values at the centroids along the stencil
                fs_d_ctd=f_old(:,s_cell_1(:,2));
                fs_u_ctd=f_old(:,s_cell_1(:,3));
            elseif face_type==3 || face_type==4
                %% Create the pdf values at the stencil points
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                %% Create the pdf values at the centroids along the stencil
                fs_d_ctd=f_old(:,s_cell_234(:,2));
                fs_u_ctd=f_old(:,s_cell_234(:,3));
            else
                error('Wrong flag for interior faces!');
            end
            fs_dd=fs_d; % Dummy
            fs_uu=fs_u; % Dummy
            fs_dd_ctd=fs_d_ctd; % Dummy
            fs_uu_ctd=fs_u_ctd; % Dummy
        else % Other flux calculation schemes that requires all four stencil points
            if face_type==1
                %% Create the pdf values at the stencil points
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),gf,fmp,V);
                %% Create the pdf values at the centroids along the stencil
                fs_d_ctd=f_old(:,s_cell_1(:,2));
                fs_u_ctd=f_old(:,s_cell_1(:,3));
                if s_cell_1(:,1)==0
                    %% Create the pdf values at the stencil points
                    fs_intcpt=(1-s_b_intcp_1(3,1))*f_nd(:,s_b_intcp_1(1,1))+s_b_intcp_1(3,1)*f_nd(:,s_b_intcp_1(2,1));
%                     ND1=NODE{s_b_intcp_1(1,1)};
%                     ND2=NODE{s_b_intcp_1(2,1)};
%                     ratio=s_b_intcp_1(3,1);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_dd=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_1(:,2)),dis(Coor_intcp,s_stcl_coord_1(:,3)),fs_intcpt,fs_d,fs_u,-dis(Coor_intcp,s_stcl_coord_1(:,1)),fe);
                    fs_dd=2*fs_intcpt-fs_d;
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=fs_dd;
                else
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),gf,fmp,V);
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_1(:,1));
                end
                if s_cell_1(:,4)==0
                    %% Create the pdf values at the stencil points
                    fs_intcpt=(1-s_b_intcp_1(3,4))*f_nd(:,s_b_intcp_1(1,4))+s_b_intcp_1(3,4)*f_nd(:,s_b_intcp_1(2,4));
%                     ND1=NODE{s_b_intcp_1(1,4)};
%                     ND2=NODE{s_b_intcp_1(2,4)};
%                     ratio=s_b_intcp_1(3,4);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_uu=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_1(:,3)),dis(Coor_intcp,s_stcl_coord_1(:,2)),fs_intcpt,fs_u,fs_d,-dis(Coor_intcp,s_stcl_coord_1(:,4)),fe);
                    fs_uu=2*fs_intcpt-fs_u;
                    %% Create the pdf values at the centroids along the stencil
                    fs_uu_ctd=fs_uu;
                else
                    %% Create the pdf values at the stencil points
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),gf,fmp,V);
                    %% Create the pdf values at the centroids along the stencil
                    fs_uu_ctd=f_old(:,s_cell_1(:,4));
                end
            elseif face_type==3 || face_type==4
                %% Create the pdf values at the stencil points
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                %% Create the pdf values at the centroids along the stencil
                fs_d_ctd=f_old(:,s_cell_234(:,2));
                fs_u_ctd=f_old(:,s_cell_234(:,3));
                if s_cell_234(:,1)==0
                    %% Create the pdf values at the stencil points
                    fs_intcpt=(1-s_b_intcp_234(3,1))*f_nd(:,s_b_intcp_234(1,1))+s_b_intcp_234(3,1)*f_nd(:,s_b_intcp_234(2,1));
%                     ND1=NODE{s_b_intcp_234(1,1)};
%                     ND2=NODE{s_b_intcp_234(2,1)};
%                     ratio=s_b_intcp_234(3,1);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_dd=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_intcpt,fs_d,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_dd=2*fs_intcpt-fs_d;
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=fs_dd;
                else
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                end
                if s_cell_234(:,4)==0
                    %% Create the pdf values at the stencil points
                    fs_intcpt=(1-s_b_intcp_234(3,4))*f_nd(:,s_b_intcp_234(1,4))+s_b_intcp_234(3,4)*f_nd(:,s_b_intcp_234(2,4));
%                     ND1=NODE{s_b_intcp_234(1,4)};
%                     ND2=NODE{s_b_intcp_234(2,4)};
%                     ratio=s_b_intcp_234(3,4);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_uu=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_intcpt,fs_u,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_uu=2*fs_intcpt-fs_u;
                    %% Create the pdf values at the centroids along the stencil
                    fs_uu_ctd=fs_uu;
                else
                    %% Create the pdf values at the stencil points
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    %% Create the pdf values at the centroids along the stencil
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                end
            else
                error('Wrong flag for interior faces!');
            end
        end
        %% calculate the PDFs at the face center
        if face_type==1 || face_type==4
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_14,latt_dir_marker,dt,face_normal_lattice,interior,fupd,ftvd);
        elseif face_type==3
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,interior,fupd,ftvd);
        else
            error('Logic error!');
        end
    elseif face_bc_flag>0 % Outer boundary face
        if face_bc_flag==1 % Periodic
            %% Calculate the PDFs at each stencil points
            if periodic<=1 % 1st-order upwind and LW
                %% Create the pdf values at the stencil points
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_d_ctd=f_old(:,s_cell_234(:,2));
                fs_u_ctd=f_old(:,s_cell_234(:,3));
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                %% Create the pdf values at the stencil points
                fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=f_old(:,s_cell_234(:,1));
                fs_d_ctd=f_old(:,s_cell_234(:,2));
                fs_u_ctd=f_old(:,s_cell_234(:,3));
                fs_uu_ctd=f_old(:,s_cell_234(:,4));
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,periodic,fupd,ftvd);
        elseif face_bc_flag==20 || face_bc_flag==21 % Velocity inlet or pressure inlet
            %% Calculate the PDFs at each stencil points
            if inlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
%                     ND1=NODE{s_b_intcp_234(1,2)};
%                     ND2=NODE{s_b_intcp_234(2,2)};
%                     ratio=s_b_intcp_234(3,2);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_d=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,4)),fs_intcpt,fs_u,fs_uu,-dis(Coor_intcp,s_stcl_coord_234(:,2)),fe);
%                     
%                     fs_dd=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,2)),0,dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_d,fs_intcpt,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_d=2*fs_intcpt-fs_u;
%                     fs_dd=2*fs_d-fs_intcpt;
                    fs_dd=fs_d;
%                     fs_dd=fs_d; % A even better scheme
%                     fs_dd=fs_intcpt;
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_d;
                    fs_dd_ctd=fs_dd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
%                     ND1=NODE{s_b_intcp_234(1,3)};
%                     ND2=NODE{s_b_intcp_234(2,3)};
%                     ratio=s_b_intcp_234(3,3);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_u=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,1)),fs_intcpt,fs_d,fs_dd,-dis(Coor_intcp,s_stcl_coord_234(:,3)),fe);
%                     
%                     fs_uu=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,3)),0,dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_u,fs_intcpt,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_u=2*fs_intcpt-fs_d;
%                     fs_uu=2*fs_u-fs_intcpt;
                    fs_uu=fs_u;         
%                     fs_uu=fs_u; % A better scheme
%                     fs_uu=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                    fs_uu_ctd=fs_uu;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,inlet,fupd,ftvd);
        elseif face_bc_flag==30 || face_bc_flag==31 % Velocity outlet or pressure outlet
            %% Calculate the PDFs at each stencil points
            if outlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy                
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
%                     ND1=NODE{s_b_intcp_234(1,2)};
%                     ND2=NODE{s_b_intcp_234(2,2)};
%                     ratio=s_b_intcp_234(3,2);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_d=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,4)),fs_intcpt,fs_u,fs_uu,-dis(Coor_intcp,s_stcl_coord_234(:,2)),fe);
%                     
%                     fs_dd=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,2)),0,dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_d,fs_intcpt,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_d=2*fs_intcpt-fs_u;
%                     fs_dd=2*fs_d-fs_intcpt;
                    fs_dd=fs_d; % A even better scheme
%                     fs_dd=fs_intcpt;
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_d;
                    fs_dd_ctd=fs_dd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
%                     ND1=NODE{s_b_intcp_234(1,3)};
%                     ND2=NODE{s_b_intcp_234(2,3)};
%                     ratio=s_b_intcp_234(3,3);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_u=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,1)),fs_intcpt,fs_d,fs_dd,-dis(Coor_intcp,s_stcl_coord_234(:,3)),fe);
%                     
%                     fs_uu=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,3)),0,dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_u,fs_intcpt,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_u=2*fs_intcpt-fs_d;
%                     fs_uu=2*fs_u-fs_intcpt;
                    fs_uu=fs_u; % A better scheme
%                     fs_uu=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                    fs_uu_ctd=fs_uu;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,outlet,fupd,ftvd);
        elseif face_bc_flag==4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
%                     ND1=NODE{s_b_intcp_234(1,2)};
%                     ND2=NODE{s_b_intcp_234(2,2)};
%                     ratio=s_b_intcp_234(3,2);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_d=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,4)),fs_intcpt,fs_u,fs_uu,-dis(Coor_intcp,s_stcl_coord_234(:,2)),fe);
%                     
%                     fs_dd=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,2)),0,dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_d,fs_intcpt,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_d=2*fs_intcpt-fs_u;
%                     fs_dd=2*fs_d-fs_intcpt;
                    fs_dd=fs_d; % A even better scheme
%                     fs_dd=fs_intcpt;
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_d;
                    fs_dd_ctd=fs_dd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
%                     ND1=NODE{s_b_intcp_234(1,3)};
%                     ND2=NODE{s_b_intcp_234(2,3)};
%                     ratio=s_b_intcp_234(3,3);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_u=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,1)),fs_intcpt,fs_d,fs_dd,-dis(Coor_intcp,s_stcl_coord_234(:,3)),fe);
%                     
%                     fs_uu=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,3)),0,dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_u,fs_intcpt,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_u=2*fs_intcpt-fs_d;
%                     fs_uu=2*fs_u-fs_intcpt;
                    fs_uu=fs_u; % A better scheme
%                     fs_uu=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                    fs_uu_ctd=fs_uu;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,s_wall,fupd,ftvd);
        elseif face_bc_flag==5
            %% Calculate the PDFs at each stencil points
            if m_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
%                     ND1=NODE{s_b_intcp_234(1,2)};
%                     ND2=NODE{s_b_intcp_234(2,2)};
%                     ratio=s_b_intcp_234(3,2);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_d=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,4)),fs_intcpt,fs_u,fs_uu,-dis(Coor_intcp,s_stcl_coord_234(:,2)),fe);
%                     
%                     fs_dd=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,2)),0,dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_d,fs_intcpt,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_d=2*fs_intcpt-fs_u;
%                     fs_dd=2*fs_d-fs_intcpt;
                    fs_dd=fs_d; % A better scheme
%                     fs_dd=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_d;
                    fs_dd_ctd=fs_dd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
%                     ND1=NODE{s_b_intcp_234(1,3)};
%                     ND2=NODE{s_b_intcp_234(2,3)};
%                     ratio=s_b_intcp_234(3,3);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_u=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,1)),fs_intcpt,fs_d,fs_dd,-dis(Coor_intcp,s_stcl_coord_234(:,3)),fe);
%                     
%                     fs_uu=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,3)),0,dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_u,fs_intcpt,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_u=2*fs_intcpt-fs_d;
%                     fs_uu=2*fs_u-fs_intcpt;
                    fs_uu=fs_u; % A better scheme
%                     fs_uu=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                    fs_uu_ctd=fs_uu;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,m_wall,fupd,ftvd);
        elseif face_bc_flag==6
            %% Calculate the PDFs at each stencil points
            if well_dp<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                    fs_dd=fs_d; % Zero gradient
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_u_ctd;
                    fs_dd_ctd=fs_d_ctd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                    fs_uu=fs_u; % Zero gradient
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_d_ctd;
                    fs_uu_ctd=fs_u_ctd;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,well_dp,fupd,ftvd);
        elseif face_bc_flag==8
            %% Calculate the PDFs at each stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),gf,fmp,V);
                    
                    fs_c=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    
%                     fs_c=(fs_u*s_dis(1,4)-fs_uu*s_dis(1,3))/(s_dis(1,4)-s_dis(1,3)); % Explode
                    fs_d=2*fs_c-fs_u; % Linear profile
                    fs_dd=2*fs_d-fs_c; % Linear profile
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_234(:,3));
                    fs_uu_ctd=f_old(:,s_cell_234(:,4));
                    fs_d_ctd=fs_u_ctd;
                    fs_dd_ctd=fs_d_ctd;
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),gf,fmp,V);
                    
                    fs_c=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    
%                     fs_c=(fs_d*s_dis(1,1)-fs_dd*s_dis(1,2))/(s_dis(1,1)-s_dis(1,2)); % Explode
                    fs_u=2*fs_c-fs_d; % Linear profile
                    fs_uu=2*fs_u-fs_c; % Linear profile
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_234(:,1));
                    fs_d_ctd=f_old(:,s_cell_234(:,2));
                    fs_u_ctd=fs_d_ctd;
                    fs_uu_ctd=fs_u_ctd;
                else
                    error('Logic error!');
                end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_23,latt_dir_marker,dt,face_normal_lattice,iof,fupd,ftvd);
        else
            error('Wrong flag for face on outer boundaries!');
        end
    elseif face_bc_flag<0 % Inner boundary face
        if face_bc_flag==-4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*f_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*f_nd(:,s_b_intcp_1(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_1(:,3));
                    fs_d_ctd=fs_d;
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*f_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*f_nd(:,s_b_intcp_1(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    %% Create the pdf values at the centroids along the stencil
                    fs_d_ctd=f_old(:,s_cell_1(:,2));
                    fs_u_ctd=fs_u;
                else
                    error('Logic error!');
                end
                %% Create the pdf values at the stencil points
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
                %% Create the pdf values at the centroids along the stencil
                fs_dd_ctd=fs_d_ctd; % Dummy
                fs_uu_ctd=fs_u_ctd; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    %% Create the pdf values at the stencil points
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),gf,fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*f_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*f_nd(:,s_b_intcp_1(2,2));
%                     ND1=NODE{s_b_intcp_234(1,2)};
%                     ND2=NODE{s_b_intcp_234(2,2)};
%                     ratio=s_b_intcp_234(3,2);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_d=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,3)),dis(Coor_intcp,s_stcl_coord_234(:,4)),fs_intcpt,fs_u,fs_uu,-dis(Coor_intcp,s_stcl_coord_234(:,2)),fe);
%                     
%                     fs_dd=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,2)),0,dis(Coor_intcp,s_stcl_coord_234(:,3)),fs_d,fs_intcpt,fs_u,-dis(Coor_intcp,s_stcl_coord_234(:,1)),fe);
                    fs_d=2*fs_intcpt-fs_u;
%                     fs_dd=2*fs_d-fs_intcpt;
                    fs_dd=fs_d; % A better scheme
%                     fs_dd=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_u_ctd=f_old(:,s_cell_1(:,3));
                    fs_uu_ctd=f_old(:,s_cell_1(:,4));
                    fs_d_ctd=fs_d;
                    fs_dd_ctd=fs_dd;
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    %% Create the pdf values at the stencil points
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),gf,fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),gf,fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*f_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*f_nd(:,s_b_intcp_1(2,3));
%                     ND1=NODE{s_b_intcp_234(1,3)};
%                     ND2=NODE{s_b_intcp_234(2,3)};
%                     ratio=s_b_intcp_234(3,3);
%                     Coor_intcp=ND1{3}+ratio*(ND2{3}-ND1{3});
%                     fs_u=inter_extrp(0,dis(Coor_intcp,s_stcl_coord_234(:,2)),dis(Coor_intcp,s_stcl_coord_234(:,1)),fs_intcpt,fs_d,fs_dd,-dis(Coor_intcp,s_stcl_coord_234(:,3)),fe);
%                     
%                     fs_uu=inter_extrp(-dis(Coor_intcp,s_stcl_coord_234(:,3)),0,dis(Coor_intcp,s_stcl_coord_234(:,2)),fs_u,fs_intcpt,fs_d,-dis(Coor_intcp,s_stcl_coord_234(:,4)),fe);
                    fs_u=2*fs_intcpt-fs_d;
%                     fs_uu=2*fs_u-fs_intcpt;
                    fs_uu=fs_u; % A better scheme
%                     fs_uu=fs_intcpt; % A even better scheme
                    %% Create the pdf values at the centroids along the stencil
                    fs_dd_ctd=f_old(:,s_cell_1(:,1));
                    fs_d_ctd=f_old(:,s_cell_1(:,2));
                    fs_u_ctd=fs_u;
                    fs_uu_ctd=fs_uu;
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,ls_14,latt_dir_marker,dt,face_normal_lattice,s_wall,fupd,ftvd);
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