FInr=9; % Flux flag for interior face
FPdc=9; % Flux flag for periodic outer boundary face
FInt=9; % Flux flag for intlet face
FOut=9; % Flux flag for outlet face
FStw=9; % Flux flag for stationary wall face
FMow=9; % Flux flag for moving wall face
FWdp=9; % Flux flag for well developed face


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


FInr=1; % Flux flag for interior face
FPdc=1; % Flux flag for periodic outer boundary face
FInt=1; % Flux flag for intlet face
FOut=1; % Flux flag for outlet face
FStw=1; % Flux flag for stationary wall face
FMow=1; % Flux flag for moving wall face
FWdp=1; % Flux flag for well developed face


top='Periodic';
bottom='Periodic';
%%%% Boundary condition
[CELL,NODE,FACE,FPDC]=bc(CELL,M,NODE,N,FACE,O,N_I,N_I_N,N_L,N_H,X1,X2,Y1,Y2,h,dx,dy,FM,V1,V2,top,right,bottom,left,immersed);