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
Tau=0.006;
Tau_t=0.003;
%% Define velocity decomposition of the lattice
[V1,Mew1]=vlatt(d,qh,Tau);
[V2,Mew2]=vlatt(d,qt,Tau_t);
%% Defining lattice weighting factor for equilibrium pdf
[wh,wt]=weight(d,qh,qt);
V=V1;
%% Create stencil infomation
[CELL,NODE,FACE]=stencil(CELL,NODE,FACE,V1,V2,N_I,N_I_N,N_L,N_H,h,dx,dy,X,Y,X1,X2,Y1,Y2,FM);
%% Define physical boundary conditions % 1% computation
top='Stationary Wall';
right='Periodic';
bottom='Stationary Wall';
left='Periodic';
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
U_m=[0.01/sqrt(3),0,0.01,0;0,0.1,0,0.2]; % Moving wall velocity [top,right,bottom,left]

U_in=[0,0.01,0.01,0.2/3;2,0,0,0]; % Inlet velocity [top,right,bottom,left]
F_u_in=[0,0,0,0;0,0,0,0]; % Variable distribution type for inlet velocity [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic

U_out=[0,0.01,0.01,0.01;2,0,0,0]; % Outlet velocity [top,right,bottom,left]
F_u_out=[2,0,0,0;2,0,0,0]; % Variable distribution type for outlet velocity [top,right,bottom,left]; 0 --- constant; 1 --- Linear; 2 --- Parabolic
% Thermal
T_bc=[1/3,1/3,1/3,1/3,2/3];  % Temperature [top,right,bottom,left,immersed]
TG_bc=[1,0,0,0,0;0,0,0,0,0]; % Temperature gradient [top,right,bottom,left,immersed]
%
FCO=[0,0,0,0]; %%%% Flag for variables at corner nodes. [top right,bottom right,bottom left,top left]
% FCO=0----Succeed the smaller value from neighbor boundary node;
% FCO=1----Succeed the larger value from neighbor boundary node;
% FCO=2----0.5*(max+min);
%%%%Write Macro Variables into boundary nodes
NODE=macro_bc(X1,X2,Y1,Y2,N_I,N_H,N_L,NODE,Rho_ref,Rho_in,F_rho_in,Rho_out,F_rho_out,U_r,U_m,U_in,F_u_in,U_out,F_u_out,T_bc,TG_bc,FCO);
%% Control Flags
%%%% Solver-related flags
FF=0; %%%% Flag for Flow type. FF=0----Normal channel flow; FF=1----Taylor vortex flow; FF=2----Perturbation flow
FFC=1; %%%% Flag for force term. FFC=0---There is no force term; FFC=1---force term
FS=0;  %%%% Splitting flag, 1----Collision first; 2----Advection first; 3----Strang splitting with advection in the middle; 4----Strang splitting with collision in the middle;
FT=1; %%%% Time-matching scheme. FT=0----Euler; FT=1----Implicit Euler; FT=2----AB;FT=3----Runge Kutta
FD=0; %%%% Flag of density used for equilibrium and macro. FD=0----Local density; FD=1----Reference density
FUPD=0; %% Flag for the flag for wether use upwind scheme to calculate the flux
FTVD=0; %%% Flag for TVD scheme
FDC=0; % Flag for decaying of f_neq in pdf_bc_h
FMP=1; % Flag for mapping, FMP=1---First-order; FMP=2---Second-order
FE=0; % FE is the flag for extrapolation of pdf along the stencil for the ghost stencil points on boundary. FE=0---linear extrapolation;FE=1---cubic extrapolation
FTHM=0; % Flag for different thermal model. FTHM=0---Passive Scalar model; FTHM=1---Double Distribution Model; FTHM=2---High Order Lattice model
FHYD_solver_switch=0; % The switch to turn on and off hydrodynamic solver 0---off; 1---on
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
FPdc=1; % Flux flag for periodic outer boundary face
FInt=1; % Flux flag for intlet face
FOut=1; % Flux flag for outlet face
FStw=1; % Flux flag for stationary wall face
FMow=1; % Flux flag for moving wall face
FWdp=1; % Flux flag for well developed face
FIof=1; % Flux flag for well developed face
%% Other variables
tt=1; % Time step counter
TT(1)=tt; % The initialization of history time step counter
W=1; %%%% Weighting factor for triangle
wl=0; %%%% Weighting factor for flux
wlb=0; %%%% Weighting factor for flux on boundary
wc=0.5; %%%% Weighting factor for extrapolation of corner nodes
Body_force=[0;-0.04];
T_ref=1/3;
%%%% Determing the size of dt
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
dt2=0.4*Tau;
dt=min(dt1,dt2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%PDF Initialization%%%%%%%%%%%%%%%%%%%%%%%%%
ini
%%
%%%%Iteration starts
solver
%%
plt