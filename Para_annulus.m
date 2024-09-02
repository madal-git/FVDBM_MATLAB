feature accel on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Parameters Setup%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FM==0 % FM----FM=0----IRT mesh; FM=1----Random mesh
    N_I_N=0; % The IRT mesh has no immersed boundary nodes.
    dx=0;
    dy=0;
else
    h=0;
end
FTH=0; %%%% Flag for thermal model. 0---Isothermal; 1---Thermal
FCD=1; %%% Flag for whether the problem is purely heat conduction. 0----No; 1----Yes
d=2;
qh=9;
qt=9;
if d~=2
    error('Sorry, the current  model could only do 2-D simulations!');
end
if qh==37
    FTH=1; % Flag for thermal model
end

if qh~=7 && qh~=9 && qh~=13 && qh~=37
    error('Sorry, the current model could only handle D2Q7, D2Q9, D2Q13 or D2Q37 lattice for hydrodynamic pdf!');
end
% Tau=sqrt(3)/500;
Tau=0.006;
Tau_t=0.003;
%% Define velocity decomposition of the lattice
[V1,Mew1]=vlatt(d,qh,Tau);
[V2,Mew2]=vlatt(d,qt,Tau_t);
%% Defining lattice weighting factor for equilibrium pdf
[wh,wt]=weight(d,qh,qt);
V=V1;
%% Determing the size of dt
Size_min=(X2-X1)*(Y2-Y1);
for r=1:M
    P=CELL{r};
    if Size_min<P{6}
        ;
    else
        Size_min=P{6};
    end
end
h_min=sqrt(Size_min*2);
dt1=0.18*h_min;
dt2=4*Tau;
dt=min(dt1,dt2);
% dt=0.1*Tau;
%% Create stencil infomation
[CELL,NODE,FACE]=stencil(CELL,NODE,FACE,V1,V2,N_I,N_I_N,N_L,N_H,h,dx,dy,X,Y,X1,X2,Y1,Y2,dt,FM);
%% Define physical boundary conditions % 1% computation
top='Stationary Wall';
right='Stationary Wall';
bottom='Stationary Wall';
left='Stationary Wall';
immersed='Stationary Wall';
%%%% Boundary condition
[CELL,NODE,FACE,FPDC]=bc(CELL,M,NODE,N,FACE,O,N_I,N_I_N,N_L,N_H,X1,X2,Y1,Y2,h,dx,dy,FM,V1,V2,top,right,bottom,left,immersed);
%% Boundary Macro Variables Definition %%%%%%%%%
% Hydro
Rho_ref=2; % Reference density

Rho_in=[2,2,2,2]; % Inlet density [top,right,bottom,left]
F_rho_in=[0,0,0,0]; % Variable distribution type for inlet density [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic

Rho_out=[2,1.95,2,2]; % Outlet density [top,right,bottom,left]
F_rho_out=[2,0,0,0]; % Variable distribution type for outlet density [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic

U_r=[0,0,0,0;0,0,0,0]; % Reference velocity
% U_m=[0.1/sqrt(3),0,0,0;0,0,0,0]; % Moving wall velocity [top,right,bottom,left]
U_m=[0.1,0,0,0;0,0,0,0]; % Moving wall velocity [top,right,bottom,left]


U_in=[0,0.01,0.01,0.2/3;2,0,0,0]; % Inlet velocity [top,right,bottom,left]
F_u_in=[0,0,0,0;0,0,0,0]; % Variable distribution type for inlet velocity [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic

U_out=[0,0.01,0.01,0.01;2,0,0,0]; % Outlet velocity [top,right,bottom,left]
F_u_out=[2,0,0,0;2,0,0,0]; % Variable distribution type for outlet velocity [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic
% Thermal
T_bc=[1/3,1/3,1/3,1/3,2/3];  % Temperature [top,right,bottom,left,immersed]
TG_bc=[1,0,0,0,0;0,0,0,0,0]; % Temperature gradient [top,right,bottom,left,immersed]
F_T=[0,0,0,0;0,0,0,0]; % Variable distribution type for temperature [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic
F_TG=[0,0,0,0;0,0,0,0]; % Variable distribution type for temperature gradient [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic
%
FCO=[0,0,0,0]; %%%% Flag for variables at corner nodes. [top right,bottom right,bottom left,top left]
% FCO=0----Succeed the smaller value from neighbor boundary node;
% FCO=1----Succeed the larger value from neighbor boundary node;
% FCO=2----0.5*(max+min);
%%%%Write Macro Variables into boundary nodes
NODE=macro_bc(X1,X2,Y1,Y2,N_I,N_H,N_L,NODE,Rho_ref,Rho_in,F_rho_in,Rho_out,F_rho_out,U_r,U_m,U_in,F_u_in,U_out,F_u_out,T_bc,TG_bc,F_T,F_TG,FCO);
%% Control Flags
%%%% Solver-related flags
FF=0; %%%% Flag for Flow type. FF=0----Normal channel flow; FF=1----Taylor vortex flow; FF=2----Perturbation flow
FFC=0; %%%% Flag for force term. FFC=0---There is no force term; FFC=1---force term
FS=0;  %%%% Splitting flag, 1----Collision first; 2----Advection first; 3----Strang splitting with advection in the middle; 4----Strang splitting with collision in the middle;
FT=6; %%%% Time-matching scheme. FT=0----Euler; FT=1----Implicit Euler; FT=2----AB;FT=3----Runge Kutta
FD=0; %%%% Flag of density used for equilibrium and macro. FD=0----Local density; FD=1----Reference density
FUPD=0; %% Flag for the flag for wether use upwind scheme to calculate the flux
FTVD=0; %%% Flag for TVD scheme
FDC=0; % Flag for decaying of f_neq in pdf_bc_h
FMP=1; % Flag for mapping, FMP=1---First-order; FMP=2---Second-order
FE=0; % FE is the flag for extrapolation of pdf along the stencil for the ghost stencil points on boundary. FE=0---linear extrapolation;FE=1---cubic extrapolation
FTHM=0; % Flag for different thermal model. FTHM=0---Passive Scalar model; FTHM=1---Double Distribution Model; FTHM=2---High Order Lattice model
FHYD_solver_switch=1; % The switch to turn on and off hydrodynamic solver 0---off; 1---on
FMRT=0; % 0---SRT model; 1---TRT model; 2---MRT model
if FMRT==0
    ;
elseif FMRT==1
    error('Temporarily unavailable!');
elseif FMRT==2
    [M_h,M_h_inv] = mrt_2d_trans_matrix (d,qh);
    S_h = mrt_2d_relax_matrix (d,qh,Tau);
else
    error('Wrong flag for collision model!');
end
%%%% Define SLIC
FSLIC=0; % This is the flag for SLIC scheme.0---1st-order interpolation; 1---2nd-order interpolation
if FSLIC==1
    [CELL,NODE,FACE]=in_cell_mapping_SLIC_create_or_update(CELL,NODE,FACE,X1,X2,Y1,Y2,V,dt);
end
FTRC=0; % 0--- Turn off tracer plot; 1--- Turn on tracer plot
if FTRC==1
    Tracer_plot_interval=800;
    Num_tracer_x=15;
    Num_tracer_y=15;
end
FANI=0; % Flag for other animation plot
if FANI==1
    Anim_plot_interval=800;
end
%%%%Define scheme flag for boundary nodal pdf update (nodal boundary scheme)
NPdc=131;
NInt=2010; % 2010 for U inlet; 2100 for Rho inlet; 2099 for U and Rho Inlet
NOut=3100; % 3100 for extrapolated u and v; 3199 for extrapolated u, v=0
NStw=410;
NMow=510;
NWdp=620;
NIof=800;
NC70=7000;
NC71=7110;
NC72=7210;
NC73=7310; % 7310 for only velocity presribed; 7399 for both density and velocity prescribed
NC74=7400;
NC75=7500;
NC76=7620;
NC77=7700;
NC78=7800;
NC79=7900;
NC80=8000;
for i=1:N
    ND=NODE{i};
    if ND{2}==20
        if NInt==2099
            Nd_bc=ND{21};
            Nd_bc(1,1)=Rho_in(4);
            ND{21}=Nd_bc;
        end
    end
    if ND{2}==73
        if NC73==7399
            Nd_bc=ND{21};
            Nd_bc(1,1)=Rho_in(4);
            ND{21}=Nd_bc;
        end
    end
    NODE{i}=ND;
end
%%%%Define scheme flag for flux calculation
% The following flags have the same rule: 
% 0---1st-order upwind (FOU)
% 1---Lax-Wendroff (LW)
% 2---2nd-order upwind (SOU)
% 3---TVD
% 4---QUICK
% 5---QUICKEST
FInr=1; % Flux flag for interior face
FPdc=FInr; % Flux flag for periodic outer boundary face
FInt=FInr; % Flux flag for intlet face
FOut=FInr; % Flux flag for outlet face
FStw=FInr; % Flux flag for stationary wall face
FMow=FInr; % Flux flag for moving wall face
FWdp=FInr; % Flux flag for well developed face
FIof=FInr; % Flux flag for well developed face
%% Other variables
tt=1; % Time step counter
TT(1)=tt; % The initialization of history time step counter
W=1; %%%% Weighting factor for triangle
wl=0; %%%% Weighting factor for flux
wlb=0; %%%% Weighting factor for flux on boundary
wc=0.5; %%%% Weighting factor for extrapolation of corner nodes
Body_force=[0;-0.04];
T_ref=1/3;

Vol=zeros(qh,M);
for r=1:M
    P=CELL{r};
    Vol(:,r)=P{6};
end
t_stop=100;

%% Make all changes due to the created outer boundary
%% Determine the inner and outer boundary radius
plt_mesh;
Outer_boundary_node=0;
Radius_outer_boundary_node=0;
Outer_boundary_node_counter=0;
hold on
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if dis(coord,[1;1])>0.66 && dis(coord,[1;1])<0.68
        Outer_boundary_node_counter=Outer_boundary_node_counter+1;
        Outer_boundary_node(Outer_boundary_node_counter)=i;
        Radius_outer_boundary_node(Outer_boundary_node_counter)=dis(coord,[1;1]);
        plot(coord(1),coord(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        hold on
    end
end
Radius_outer=sum(Radius_outer_boundary_node)/Outer_boundary_node_counter;

Inner_boundary_node=0;
Radius_inner_boundary_node=0;
Inner_boundary_node_counter=0;
hold on
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if dis(coord,[1;1])>0.32 && dis(coord,[1;1])<0.34
        Inner_boundary_node_counter=Inner_boundary_node_counter+1;
        Inner_boundary_node(Inner_boundary_node_counter)=i;
        Radius_inner_boundary_node(Inner_boundary_node_counter)=dis(coord,[1;1]);
        plot(coord(1),coord(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
        hold on
    end
end
Radius_inner=sum(Radius_inner_boundary_node)/Inner_boundary_node_counter;

%% NODE
for i=1:Outer_boundary_node_counter
    ND=NODE{Outer_boundary_node(i)};
    ND{2}=5; % Moving wall
    NUM_STAR_CELL=ND{4};
    if NUM_STAR_CELL~=6
        error('Logic Error-NUM_STAR_CELL~=6!');
    end
    ORDER_STAR_CELL=ND{5};
    DIS_STAR_CELL=ND{6};
    %% cut the half that is outside of outer boundary
    NUM_STAR_CELL_new=0;
    ORDER_STAR_CELL_new=0;
    DIS_STAR_CELL_new=0;
    DISRE_STAR_CELL_new=0;
    for r=1:NUM_STAR_CELL
        CL=CELL{ORDER_STAR_CELL(r)};
        if dis(CL{5},[1;1])<Radius_outer
            NUM_STAR_CELL_new=NUM_STAR_CELL_new+1;
            ORDER_STAR_CELL_new(NUM_STAR_CELL_new)=ORDER_STAR_CELL(r);
            DIS_STAR_CELL_new(NUM_STAR_CELL_new)=DIS_STAR_CELL(r);
            DISRE_STAR_CELL_new=DISRE_STAR_CELL_new+1/DIS_STAR_CELL(r);
        end
    end
    if NUM_STAR_CELL_new~=3
        error('Logic Error-NUM_STAR_CELL_new~=3!');
    end
    ND{4}=NUM_STAR_CELL_new;
    ND{5}=ORDER_STAR_CELL_new;
    ND{6}=DIS_STAR_CELL_new;
    ND{7}=DISRE_STAR_CELL_new;
    NODE{Outer_boundary_node(i)}=ND;
end
% Check
% for i=1:Outer_boundary_node_counter
%     ND=NODE{Outer_boundary_node(i)};
%     if ND{2}~=5
%         error('Wrong flag for outer cylinder boundary!');
%     end
%     if ND{4}~=3
%         error('Logic Error!');
%     end
%     ND{5}
%     ND{6}
%     ND{7}
% end

%% FACE
Outer_boundary_face=0;
Outer_boundary_face_counter=0;
for i=1:O
    FC=FACE{i};
    nd1=FC{8};
    nd2=FC{9};
    if length(union(nd1,Outer_boundary_node))==Outer_boundary_node_counter && length(union(nd2,Outer_boundary_node))==Outer_boundary_node_counter
        Outer_boundary_face_counter=Outer_boundary_face_counter+1;
        Outer_boundary_face(Outer_boundary_face_counter)=i;
        Coord_nd1=FC{5};
        Coord_nd2=FC{6};
        plot([Coord_nd1(1),Coord_nd2(1)],[Coord_nd1(2),Coord_nd2(2)],'red','linewidth',1.5);
        hold on
    end
end
if Outer_boundary_face_counter~=Outer_boundary_node_counter
    error('Logic Error!');
end
%fc2,fc10,fc11,fc12,fc13,fc16,fc17,fc18,fc22, fc24
for i=1:Outer_boundary_face_counter
    FC=FACE{Outer_boundary_face(i)};
    ND1=NODE{FC{8}};
    ND2=NODE{FC{9}};
    if ND1{2}~=5 || ND2{2}~=5
        error('Logic Error!');
    else
        %fc2,fc10,fc11,fc18,fc24
        FC{2}=5; % Moving Wall
        FC{10}=5;
        FC{11}=5;
        %fc12,fc13,fc16,fc17,fc18
        neigh_upwind=FC{12};
        neigh_downwind=FC{13};
        face_stencil=FC{16};
        face_stencil_point_coord=FC{17};
        face_boundary_intercept=FC{18};
        cell_neigh_upwind=CELL{neigh_upwind(1)};
        cell_neigh_downwind=CELL{neigh_downwind(1)};
        if dis(cell_neigh_upwind{5},[1;1])<Radius_outer && dis(cell_neigh_downwind{5},[1;1])>Radius_outer
            neigh_downwind=[0,0];
            FC{13}=neigh_downwind;
            
            %fc16
            if face_stencil(1)==0 || face_stencil(2)==0
                error('Logic Error!');
            else
                face_stencil(1)=0;
                face_stencil(2)=0;
            end
            
            %fc17
            C_b=FC{7}; % The base point coordinates of the stencil
            face_stencil_point_coord(:,2)=C_b+FC{4}'*dis(C_b,face_stencil_point_coord(:,3));
            face_stencil_point_coord(:,1)=face_stencil_point_coord(:,2)+FC{4}'*dis(C_b,face_stencil_point_coord(:,2));
            
            %fc18
            if sum(face_boundary_intercept(:,1))~=0 || sum(face_boundary_intercept(:,2))~=0
                error('Logic Error!');
            else
                face_boundary_intercept(1,1)=FC{8};
                face_boundary_intercept(2,1)=FC{9};
                face_boundary_intercept(3,1)=0.5;
                
                face_boundary_intercept(1,2)=FC{8};
                face_boundary_intercept(2,2)=FC{9};
                face_boundary_intercept(3,2)=0.5;
            end
        elseif dis(cell_neigh_upwind{5},[1;1])>Radius_outer && dis(cell_neigh_downwind{5},[1;1])<Radius_outer
            neigh_upwind=[0,0];
            FC{12}=neigh_upwind;
            
            %fc16
            if face_stencil(3)==0 || face_stencil(4)==0
                error('Logic Error!');
            else
                face_stencil(3)=0;
                face_stencil(4)=0;
            end
            
            %fc17
            C_b=FC{7}; % The base point coordinates of the stencil
            face_stencil_point_coord(:,3)=C_b-FC{4}'*dis(C_b,face_stencil_point_coord(:,2));
            face_stencil_point_coord(:,4)=face_stencil_point_coord(:,3)-FC{4}'*dis(C_b,face_stencil_point_coord(:,3));
            
            %fc18
            if sum(face_boundary_intercept(:,3))~=0 || sum(face_boundary_intercept(:,4))~=0
                error('Logic Error!');
            else
                face_boundary_intercept(1,3)=FC{8};
                face_boundary_intercept(2,3)=FC{9};
                face_boundary_intercept(3,3)=0.5;
                
                face_boundary_intercept(1,4)=FC{8};
                face_boundary_intercept(2,4)=FC{9};
                face_boundary_intercept(3,4)=0.5;
            end
        else
            error('Logic Error!');
        end
        FC{16}=face_stencil;
        FC{17}=face_stencil_point_coord;
        FC{18}=face_boundary_intercept;
        
        %fc22
        fc22=zeros(1,length(FC{16}));
        fc17=FC{17};
        C_b=FC{7};
        for j=1:length(FC{16})
            fc22(1,j)=dis(C_b,fc17(:,j));
        end
        FC{22}=fc22;
        
        %fc24
        if FC{24}~=1
            error('Logic Error!');
        else
            FC{24}=2;
        end
    end
    FACE{Outer_boundary_face(i)}=FC;
end
% fc46
for i=1:Outer_boundary_face_counter
    FC=FACE{Outer_boundary_face(i)};
    
    fc16=FC{16};
    fc17=FC{17};
    fc22=FC{22};
    q1=length(V1);
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
    
    %% S{30}
    s30=zeros(2,q1);
    for k=1:q1
        %s30(:,k)=FC{7}+(-V1(:,k)*dt);
        s30(:,k)=FC{7}+(-(FC{4}'.*(FC{4}*V1(:,k)))*dt);
        % Check if in correct triangle
%         if k~=1
%             if s2(1,k)~=0
%                 cell_in=CELL{s2(1,k)};
%                 if ~in_triangle(s30(:,k),cell_in{13},cell_in{14},cell_in{15})
%                     if FC{2}==0
%                         cell_in=CELL{s2(1,k)};
%                         plot(cell_in{22},cell_in{23},cell_in{24},cell_in{25},cell_in{26},cell_in{27},'black','linewidth', 1);
%                         hold on
%                         face_center=FC{7};
%                         Advection_origin_point=s30(:,k);
%                         plot(face_center(1,1),face_center(2,1),'Marker', 'o','Markersize',4, 'MarkerFaceColor', 'red','MarkerEdgeColor', 'red');
%                         hold on
%                         plot(Advection_origin_point(1,1),Advection_origin_point(2,1),'Marker', '*','Markersize',4, 'MarkerFaceColor', 'blue','MarkerEdgeColor', 'blue');
%                         hold on
%                         error('Please decrease the size of dt!');
%                     end
%                 end
%             else % The found coordinate is out of the computational domain
%                 if FC{24}==1
%                     if s30(1,k)<=X1 || s30(1,k)>=X2 || s30(2,k)<=Y1 || s30(2,k)>=Y2
%                         error('Logic Error!');
%                     end
%                 elseif FC{24}==2
%                     if ~(s30(1,k)<=X1 || s30(1,k)>=X2 || s30(2,k)<=Y1 || s30(2,k)>=Y2)
%                         error('Logic Error!');
%                     end
%                 else
%                     error('Logic Error!');
%                 end
%             end
%         end
    end
    S{30}=s30;
    %% S{31} and S{32}
    s31=zeros(q1,1);
    s32=zeros(q1,1);
    for k=1:q1
        if L_V_1(k)>=0
            s31(k,1)=dis(s30(:,k),fc17(:,2))/dis(fc17(:,3),fc17(:,2));
            s32(k,1)=dis(s30(:,k),fc17(:,3))/dis(fc17(:,3),fc17(:,2));
        else
            s31(k,1)=dis(s30(:,k),fc17(:,3))/dis(fc17(:,3),fc17(:,2));
            s32(k,1)=dis(s30(:,k),fc17(:,2))/dis(fc17(:,3),fc17(:,2));
        end
    end
    % check
    for k=1:q1
        if single(s31(k,1)+s32(k,1))~=1
            error('Logic Error-not on a line!');
        end
        if FM==0
            if single(s31(k,1))==0.5 && single(s32(k,1))==0.5
                ;
            elseif s31(k,1)<0.5 || s32(k,1)>0.5
                error('Logic Error-wrong side of the face center!');
            end
        elseif FM==1
            ;
        else
            error('Wrong flag for mesh type');
        end
    end
    S{31}=s31;
    S{32}=s32;
    %% Fill
    FC{46}=S;
    FACE{Outer_boundary_face(i)}=FC;
end
%% CELL
Active_cell=zeros(1,M);
for i=1:M
    CL=CELL{i};
    if dis(CL{5},[1;1])<Radius_outer
        Active_cell(i)=1;
    end
end
%check
for i=1:M
    CL=CELL{i};
    if Active_cell(i)
        coord=CL{5};
        plot(coord(1),coord(2),'Marker', 'o','Markersize',3, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');
    end
end
%% Write physical bc into the newly created outer boundary
epsilon=1e-4;
w_outer=0.05; % angler velocity
counter=0;
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if ND{2}==5
        if dis(coord,[1;1])-Radius_outer>epsilon
            error('Logic Error!');
        end
        counter=counter+1;
        nd_bc=ND{21};
        nd_bc(2,1)=-(coord(2)-1)*w_outer;
        nd_bc(3,1)=(coord(1)-1)*w_outer;
        ND{21}=nd_bc;
        NODE{i}=ND;
    end
end
% Check
if counter~=Outer_boundary_node_counter
    error('Logic Error!');
else
    for i=1:N
        ND=NODE{i};
        if ND{2}==5
            nd_bc=ND{21};
            if sqrt(nd_bc(2,1)^2+nd_bc(3,1)^2)-Radius_outer*w_outer>epsilon
                error('Logic error!');
            end
        end
    end
end
%% Write physical bc into the inner boundary
epsilon=1e-4;
w_inner=-0.05; % angler velocity
counter=0;
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if ND{2}==-4
        if dis(coord,[1;1])-Radius_inner>epsilon
            error('Logic Error!');
        end
        counter=counter+1;
        nd_bc=ND{21};
        nd_bc(2,1)=-(coord(2)-1)*w_inner;
        nd_bc(3,1)=(coord(1)-1)*w_inner;
        ND{21}=nd_bc;
        NODE{i}=ND;
    end
end
% Check
if counter~=Inner_boundary_node_counter
    error('Logic Error!');
else
    for i=1:N
        ND=NODE{i};
        if ND{2}==-4
            nd_bc=ND{21};
            if abs(sqrt(nd_bc(2,1)^2+nd_bc(3,1)^2)-abs(Radius_inner*w_inner))>epsilon
                error('Logic error!');
            end
        end
    end
end

%% Velocty check
U_nd_bc_a=[0;0];
U_nd_bc=[0;0];
counter=0;
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if ND{2}==5
        counter=counter+1;
        if dis(coord,[1;1])-Radius_outer>epsilon
            error('Logic Error!');
        end
        nd_bc=ND{21};
%         ND{4}
%         ND{5}
%         ND{6}
        U_nd_bc_a(:,counter)=nd_bc(2:3);
        f_nd(:,i)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,ND,CELL,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_ref,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80);
        [Rho_nd_bc,U_nd_bc(:,counter),T_nd_bc]=macro_h(f_nd(:,i),V,Rho_ref,FD);
%         U_nd_bc-U_nd_bc_a
    end
end
plot(U_nd_bc(1,:))
hold on
plot(U_nd_bc_a(1,:))
hold on
plot(U_nd_bc(2,:))
hold on
plot(U_nd_bc_a(2,:))



U_nd_bc_a=[0;0];
U_nd_bc=[0;0];
counter=0;
for i=1:N
    ND=NODE{i};
    coord=ND{3};
    if ND{2}==-4
        counter=counter+1;
        if dis(coord,[1;1])-Radius_inner>epsilon
            error('Logic Error!');
        end
        nd_bc=ND{21};
%         ND{4}
%         ND{5}
%         ND{6}
        U_nd_bc_a(:,counter)=nd_bc(2:3);
        f_nd(:,i)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,ND,CELL,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_ref,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NIof,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77,NC78,NC79,NC80);
        [Rho_nd_bc,U_nd_bc(:,counter),T_nd_bc]=macro_h(f_nd(:,i),V,Rho_ref,FD);
%         U_nd_bc-U_nd_bc_a
    end
end
plot(U_nd_bc(1,:))
hold on
plot(U_nd_bc_a(1,:))
hold on
plot(U_nd_bc(2,:))
hold on
plot(U_nd_bc_a(2,:))