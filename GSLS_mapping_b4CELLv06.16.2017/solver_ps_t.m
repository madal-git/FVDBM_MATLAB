tic;
for i=1:1000
    if FT==1 % Implicit Euler, CI/AE
        %%%%%% Update macro temperature at cell centroid
        T=macro_t(g_old,V2,Rho);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            g_eq(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                ;
            else % nodes on boundary
                g_nd(:,l)=pdf_bc_t(N_L,N_H,N_I,Tau_t,dt,NODE,NP,g_old,g_eq,g_nd,U_nd,U,Rho_nd,Rho,V2,V1,V2,Rho_r,Rho_in,Rho_out,qt,wc,wt,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end;
        end;
        %%%%% Nodal variable update
        %%%% Calculate flux integral for all triangles
        FLX_t=flux_pdf_t(CELL,M,NODE,N,FACE,O,g_nd,g_old,g_eq,Rho,Rho_r,V2,dt,wl,wlb,wt,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FPPI,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate the extrapolated f_eq_im
        g_eq_im=2*g_eq-g_eq_old;
        %%%% Calculate the current pdf at cell centroid-Implicit scheme
        for r=1:M
            P=CELL{r};
            %% 1st Algorithm, fully implicit for collision term
            g(:,r)=(g_old(:,r)-FLX_t(:,r)*dt/P{6}+g_eq_im(:,r)*dt/Tau_t)/(1+dt/Tau_t);
            %% 2nd Algorithm, half implicit for collision term
            %                 g(:,r)=((1-dt/Tau/2)*g_old(:,r)-FLX(:,r)*dt/P{6}+(f_eq_im(:,r)+f_eq(:,r))*dt/Tau/2)/(1+dt/Tau/2);
        end
        %%%%f_eq_old
        g_eq_old=g_eq;
    else
        error('The flag for time-matching scheme is unavailable or invalid for Passive-Scalar thermal model!');
    end
    %%%% Monitored variables
    R_t(tt+1)=norm(g-g_old,'fro');
    
    TR(tt+1)=0;
    for r=1:M
        TR(tt+1)=TR(tt+1)+T(1,r);
    end
    TR(tt+1)=TR(tt+1)/M;
    
    FLX_c_t(1:qt,tt+1)=0;
    for k=1:qt
        for r=1:M
            FLX_c_t(k,tt+1)=FLX_c_t(k,tt+1)+FLX_t(k,r);
        end
    end
    FLX_c_t(1:qt,tt+1)=FLX_c_t(1:qt,tt+1)/M;
    
    %%%% Update pdf at cell centroids and nodes
    g_old=g;
    g_nd_old=g_nd;
    %%%% Timer evolution
    tt=tt+1;
    TT(tt)=tt;
end
toc;
%%%%%%Nodal pdf updated and boundary conditions applied
T=macro_t(g_old,V2,Rho);
%%%% Calculate equilibrium pdf at cell centroids
for r=1:M
    g_eq(:,r)=eqm_t(V2,U(:,r),Rho(1,r),T(1,r),qt,wt,Rho_r,FD);
end
for l=1:N
    NP=NODE{l};
    if NP{2}==0 % Interior nodes
        % if W~=1
        g_nd(:,l)=node_star(NP,g_old,0);
        % end
    else % nodes on boundary
        g_nd(:,l)=pdf_bc_t(N_L,N_H,N_I,Tau_t,dt,NODE,NP,g_old,g_eq,g_nd,U_nd,U,Rho_nd,Rho,V2,V1,V2,Rho_r,Rho_in,Rho_out,qt,wc,wt,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
    end;
end;
T_nd=macro_t(g_nd,V2,Rho_nd);

T_plt=T;
for r=1:N-N_I;
    T_plt(1,M+r)=T_nd(1,r+N_I);
end
for r=1:N_I_N;
    T_plt(1,M+N-N_I+r)=T_nd(1,r);
end
Time_per_iter=toc/i