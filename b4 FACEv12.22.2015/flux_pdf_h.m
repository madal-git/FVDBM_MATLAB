function fl=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,feq,RHO,Rho_r,V,dt,wl,wlb,wh,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,ftvd,fpdc,fd,fm,fmp)
% fl=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,feq,RHO,Rho_r,V,dt,wl,wlb,wh,V1,V2,interior,periodic,inlet,outlet,s_wall,m_wall,well_dp,ftvd,fpdc,fd,fm,fmp)
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
% 1---Lax-Wendroff (LW)
% 2---2nd-order upwind (SOU)
% 3---TVD
% 4---QUICK
% 5---QUICKEST

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

for r=1:O
    %% Basic info
    FC=FACE{r};
    face_bc_flag=FC{2};
    face_type=FC{23};
    face_normal_lattice=FC{4}*V;
    %% Choose lattice-dependent info
    if lattice==1
        ls_no_pdc=FC{45};
        ls_all_pdc=FC{47};
        ls_lr_pdc=FC{49};
        ls_tb_pdc=FC{51};
        latt_dir_marker=FC{53};
    elseif lattice==2
        ls_no_pdc=FC{46};
        ls_all_pdc=FC{48};
        ls_lr_pdc=FC{50};
        ls_tb_pdc=FC{52};
        latt_dir_marker=FC{54};
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
        s_map_coord_1=zeros(6,length(FC{16}));
        
        s_cell_234=FC{16};
        s_stcl_coord_234=FC{17};
        s_b_intcp_234=FC{18};
        s_map_id_234=FC{20};
        s_map_coord_234=zeros(6,length(FC{16}));
        
        for i=1:length(s_cell_234)
            if s_cell_234(1,i)~=0
                if s_map_id_234(1,i)==0
                    Cell1=CELL{s_map_id_234(2,i)};
                    Cell2=CELL{s_map_id_234(3,i)};
                    Cell3=CELL{s_map_id_234(4,i)};
                    s_map_coord_234(1:2,i)=Cell1{5};
                    s_map_coord_234(3:4,i)=Cell2{5};
                    s_map_coord_234(5:6,i)=Cell3{5};
                elseif s_map_id_234(1,i)==1
                    Node1=NODE{s_map_id_234(2,i)};
                    Cell2=CELL{s_map_id_234(3,i)};
                    Cell3=CELL{s_map_id_234(4,i)};
                    s_map_coord_234(1:2,i)=Node1{3};
                    s_map_coord_234(3:4,i)=Cell2{5};
                    s_map_coord_234(5:6,i)=Cell3{5};
                elseif s_map_id_234(1,i)==2
                    Node1=NODE{s_map_id_234(2,i)};
                    Node2=NODE{s_map_id_234(3,i)};
                    Cell3=CELL{s_map_id_234(4,i)};
                    s_map_coord_234(1:2,i)=Node1{3};
                    s_map_coord_234(3:4,i)=Node2{3};
                    s_map_coord_234(5:6,i)=Cell3{5};
                else
                    error('Logic error!');
                end
            end
        end
        
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_no_pdc;
    elseif fpdc==1 % Left and right boundaries
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=zeros(6,length(FC{16}));
        
        s_cell_234=FC{30};
        s_stcl_coord_234=FC{31};
        s_b_intcp_234=FC{33};
        s_map_id_234=FC{35};
        s_map_coord_234=FC{36};
        
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_lr_pdc;
    elseif fpdc==2 % Top and bottom boundaries
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=zeros(6,length(FC{16}));
        
        s_cell_234=FC{37};
        s_stcl_coord_234=FC{38};
        s_b_intcp_234=FC{40};
        s_map_id_234=FC{42};
        s_map_coord_234=FC{43};
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_tb_pdc;
    elseif fpdc==3 % All periodic boundaries
        %% Select stencil info
        s_cell_1=FC{16};
        s_stcl_coord_1=FC{17};
        s_b_intcp_1=FC{18};
        s_map_id_1=FC{20};
        s_map_coord_1=zeros(6,length(FC{16}));
        
        s_cell_234=FC{24};
        s_stcl_coord_234=FC{25};
        s_b_intcp_234=zeros(1,length(FC{16}));
        s_map_id_234=FC{28};
        s_map_coord_234=FC{29};
        %% Select lattice stencil info
        ls_14=ls_no_pdc;
        ls_23=ls_all_pdc;
    else
        error('Wrong flag for periodic boundary conditions!');
    end
    
    %% Fill in coordinates in s_map_coord_1
    for i=1:length(s_cell_1)
        if s_cell_1(1,i)~=0
            if s_map_id_1(1,i)==0
                Cell1=CELL{s_map_id_1(2,i)};
                Cell2=CELL{s_map_id_1(3,i)};
                Cell3=CELL{s_map_id_1(4,i)};
                s_map_coord_1(1:2,i)=Cell1{5};
                s_map_coord_1(3:4,i)=Cell2{5};
                s_map_coord_1(5:6,i)=Cell3{5};
            elseif s_map_id_1(1,i)==1
                Node1=NODE{s_map_id_1(2,i)};
                Cell2=CELL{s_map_id_1(3,i)};
                Cell3=CELL{s_map_id_1(4,i)};
                s_map_coord_1(1:2,i)=Node1{3};
                s_map_coord_1(3:4,i)=Cell2{5};
                s_map_coord_1(5:6,i)=Cell3{5};
            elseif s_map_id_1(1,i)==2
                Node1=NODE{s_map_id_1(2,i)};
                Node2=NODE{s_map_id_1(3,i)};
                Cell3=CELL{s_map_id_1(4,i)};
                s_map_coord_1(1:2,i)=Node1{3};
                s_map_coord_1(3:4,i)=Node2{3};
                s_map_coord_1(5:6,i)=Cell3{5};
            else
                error('Logic error!');
            end
        end
    end
    ND1=NODE{FC{8}};
    ND2=NODE{FC{9}};
    f_ndl=f_nd(:,ND1{1});
    f_ndr=f_nd(:,ND2{1});
    
    
    if face_bc_flag==0 % Interior face
        %% Calculate the PDFs at each stencil points
        if interior<=1 % 1st-order upwind and LW
            if face_type==1
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
            elseif face_type==3 || face_type==4
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
            else
                error('Wrong flag for interior faces!');
            end
            fs_dd=fs_d; % Dummy
            fs_uu=fs_u; % Dummy
        else % Other flux calculation schemes that requires all four stencil points
            if face_type==1
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                if s_cell_1(:,1)==0
                    fs_intcpt=(1-s_b_intcp_1(3,1))*f_nd(:,s_b_intcp_1(1,1))+s_b_intcp_1(3,1)*f_nd(:,s_b_intcp_1(2,1));
                    fs_dd=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),fmp,V);
                end
                if s_cell_1(:,4)==0
                    fs_intcpt=(1-s_b_intcp_1(3,4))*f_nd(:,s_b_intcp_1(1,4))+s_b_intcp_1(3,4)*f_nd(:,s_b_intcp_1(2,4));
                    fs_uu=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                else
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),fmp,V);
                end
            elseif face_type==3 || face_type==4
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                if s_cell_234(:,1)==0
                    fs_intcpt=(1-s_b_intcp_234(3,1))*f_nd(:,s_b_intcp_234(1,1))+s_b_intcp_234(3,1)*f_nd(:,s_b_intcp_234(2,1));
                    fs_dd=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                end
                if s_cell_234(:,4)==0
                    fs_intcpt=(1-s_b_intcp_234(3,4))*f_nd(:,s_b_intcp_234(1,4))+s_b_intcp_234(3,4)*f_nd(:,s_b_intcp_234(2,4));
                    fs_uu=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                else
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                end
            else
                error('Wrong flag for interior faces!');
            end
        end
%         f_stencil=[fs_dd,fs_d,fs_u,fs_uu];
        %% calculate the PDFs at the face center
        if face_type==1 || face_type==4
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_14,latt_dir_marker,dt,face_normal_lattice,interior,ftvd);
        elseif face_type==3
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,interior,ftvd);
        else
            error('Logic error!');
        end
    elseif face_bc_flag>0 % Outer boundary face
        if face_bc_flag==1 % Periodic
            %% Calculate the PDFs at each stencil points
            if periodic<=1 % 1st-order upwind and LW
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
            end
%             f_stencil=[fs_dd,fs_d,fs_u,fs_uu];
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,periodic,ftvd);
        elseif face_bc_flag==20 || face_bc_flag==21 % Velocity inlet or pressure inlet
            %% Calculate the PDFs at each stencil points
            if inlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,inlet,ftvd);
        elseif face_bc_flag==30 || face_bc_flag==31 % Velocity outlet or pressure outlet
            %% Calculate the PDFs at each stencil points
            if outlet<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,outlet,ftvd);
        elseif face_bc_flag==4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
%             f_stencil=[fs_dd,fs_d,fs_u,fs_uu];
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,s_wall,ftvd);
        elseif face_bc_flag==5
            %% Calculate the PDFs at each stencil points
            if m_wall<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,2))*f_nd(:,s_b_intcp_234(1,2))+s_b_intcp_234(3,2)*f_nd(:,s_b_intcp_234(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_234(3,3))*f_nd(:,s_b_intcp_234(1,3))+s_b_intcp_234(3,3)*f_nd(:,s_b_intcp_234(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,m_wall,ftvd);
        elseif face_bc_flag==6
            %% Calculate the PDFs at each stencil points
            if well_dp<=1 % 1st-order upwind and LW
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_234(:,1)==0 && s_cell_234(:,2)==0) && (s_cell_234(:,3)~=0 && s_cell_234(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,3),s_stcl_coord_234(:,3),s_map_id_234(:,3),s_map_coord_234(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,4),s_stcl_coord_234(:,4),s_map_id_234(:,4),s_map_coord_234(:,4),fmp,V);
                    
                    fs_d=fs_u; % Zero gradient
                    fs_dd=fs_d; % Zero gradient
                elseif (s_cell_234(:,1)~=0 && s_cell_234(:,2)~=0) && (s_cell_234(:,3)==0 && s_cell_234(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,1),s_stcl_coord_234(:,1),s_map_id_234(:,1),s_map_coord_234(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_234(:,2),s_stcl_coord_234(:,2),s_map_id_234(:,2),s_map_coord_234(:,2),fmp,V);
                    
                    fs_u=fs_d; % Zero gradient
                    fs_uu=fs_u; % Zero gradient
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_23,latt_dir_marker,dt,face_normal_lattice,well_dp,ftvd);
        else
            error('Wrong flag for face on outer boundaries!');
        end
    elseif face_bc_flag<0 % Inner boundary face
        if face_bc_flag==-4
            %% Calculate the PDFs at each stencil points
            if s_wall<=1 % 1st-order upwind and LW
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*f_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*f_nd(:,s_b_intcp_1(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*f_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*f_nd(:,s_b_intcp_1(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
                fs_dd=fs_d; % Dummy
                fs_uu=fs_u; % Dummy
            else % Other flux calculation schemes that requires all four stencil points
                if (s_cell_1(:,1)==0 && s_cell_1(:,2)==0) && (s_cell_1(:,3)~=0 && s_cell_1(:,4)~=0)
                    fs_u=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,3),s_stcl_coord_1(:,3),s_map_id_1(:,3),s_map_coord_1(:,3),fmp,V);
                    fs_uu=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,4),s_stcl_coord_1(:,4),s_map_id_1(:,4),s_map_coord_1(:,4),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,2))*f_nd(:,s_b_intcp_1(1,2))+s_b_intcp_1(3,2)*f_nd(:,s_b_intcp_1(2,2));
                    fs_d=2*fs_intcpt-fs_u; % Linear extrapolation on boundary
                    
                    fs_dd=2*fs_d-fs_intcpt; % Linear extrapolation on boundary
                elseif (s_cell_1(:,1)~=0 && s_cell_1(:,2)~=0) && (s_cell_1(:,3)==0 && s_cell_1(:,4)==0)
                    fs_dd=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,1),s_stcl_coord_1(:,1),s_map_id_1(:,1),s_map_coord_1(:,1),fmp,V);
                    fs_d=in_cell_mapping(FC,CELL,NODE,f_old,f_nd,s_cell_1(:,2),s_stcl_coord_1(:,2),s_map_id_1(:,2),s_map_coord_1(:,2),fmp,V);
                    
                    fs_intcpt=(1-s_b_intcp_1(3,3))*f_nd(:,s_b_intcp_1(1,3))+s_b_intcp_1(3,3)*f_nd(:,s_b_intcp_1(2,3));
                    fs_u=2*fs_intcpt-fs_d; % Linear extrapolation on boundary
                    
                    fs_uu=2*fs_u-fs_intcpt; % Linear extrapolation on boundary
                else
                    error('Logic error!');
                end
            end
            %% calculate the PDFs at the face center
            f_face=pdf_face(fs_dd,fs_d,fs_u,fs_uu,ls_14,latt_dir_marker,dt,face_normal_lattice,s_wall,ftvd);
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