FTH=1;
qt=9;
Tau_t=0.003;
[V2,Mew2]=vlatt(d,qt,Tau_t);
FPPI=0;
%% Defining lattice weighting factor for equilibrium pdf
[wh,wt]=weight(d,qh,qt);
T_bc=[1.001,1.1,1,1.3,1.5];  % Temperature [top,right,bottom,left,immersed]
NODE=macro_bc(X1,X2,Y1,Y2,N_I,N_H,N_L,NODE,Rho_r,Rho_in,F_rho_in,Rho_out,F_rho_out,U_r,U_m,U_in,F_u_in,U_out,F_u_out,T_bc,TG_bc,FCO);
Rho=ones(1,M)*2;
Rho_nd=ones(1,N)*2;
U=zeros(2,M);
U_nd=zeros(2,N);


a=U_m(1,1)/(Y2-Y1);
b=U_m(1,1)*(1-Y2/(Y2-Y1));

for r=1:M
    P=CELL{r};
    Coor=P{5};
    U(1,r)=a*Coor(2)+b;
end

for r=1:N
    ND=NODE{r};
    Coor=ND{3};
    U_nd(1,r)=a*Coor(2)+b;
end

U_plt=U;
U_M=zeros(1,M);
for r=1:M
    P=CELL{r};
    U_M(r)=sqrt(U(1,r)^2+U(2,r)^2);
end
for r=1:N-N_I;
    U_plt(:,M+r)=U_nd(:,r+N_I);
    U_M(M+r)=sqrt(U_nd(1,r+N_I)^2+U_nd(2,r+N_I)^2);
end

figure;
Z=griddata(XXX,YYY,U_M,Xx,Yy);
contourf(Xx,Yy,Z,100);
axis equal tight;

figure;
quiver(XX,YY,U_plt(1,1:M),U(2,1:M),10);
axis equal tight  


ini_ps