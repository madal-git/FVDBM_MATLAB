FInr=5; % Flux flag for interior face
FPdc=5; % Flux flag for periodic outer boundary face
FInt=5; % Flux flag for intlet face
FOut=5; % Flux flag for outlet face
FStw=5; % Flux flag for stationary wall face
FMow=5; % Flux flag for moving wall face
FWdp=5; % Flux flag for well developed face

FInr=7; % Flux flag for interior face
FPdc=7; % Flux flag for periodic outer boundary face
FInt=7; % Flux flag for intlet face
FOut=7; % Flux flag for outlet face
FStw=7; % Flux flag for stationary wall face
FMow=7; % Flux flag for moving wall face
FWdp=7; % Flux flag for well developed face

FInr=0; % Flux flag for interior face
FPdc=0; % Flux flag for periodic outer boundary face
FInt=0; % Flux flag for intlet face
FOut=0; % Flux flag for outlet face
FStw=0; % Flux flag for stationary wall face
FMow=0; % Flux flag for moving wall face
FWdp=0; % Flux flag for well developed face

FInr=9; % Flux flag for interior face
FPdc=9; % Flux flag for periodic outer boundary face
FInt=9; % Flux flag for intlet face
FOut=9; % Flux flag for outlet face
FStw=9; % Flux flag for stationary wall face
FMow=9; % Flux flag for moving wall face
FWdp=9; % Flux flag for well developed face
FF=1
dt=dt/2
dt/Tau



U_0*sqrt(3)
dt/Tau
0.5*t_c/dt


FInr=5; % Flux flag for interior face
FPdc=5; % Flux flag for periodic outer boundary face
FInt=5; % Flux flag for intlet face
FOut=5; % Flux flag for outlet face
FStw=5; % Flux flag for stationary wall face
FMow=5; % Flux flag for moving wall face
FWdp=5; % Flux flag for well developed face
% FT=4
FF=1
C=CELL{1};
dx=sqrt(2*C{6});
dt=0.1*dx;

Tau=dt/1.5;
[V1,Mew1]=vlatt(d,qh,Tau);
[V2,Mew2]=vlatt(d,qt,Tau);
dt/Tau
dt/dx
