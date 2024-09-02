tic;
i=0;
while log10(R(end))>-6.404
    i=i+1;
    if mod(i,1000)==0
        disp(['The current run has reached ' num2str(i) ' iterations!']);
    end
    if FT==0 % Explicit Euler
        if FS==0
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    if W~=1 && FMP>1
                        %f_nd(:,l)=node_star(NP,f_old,0);
                    end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end
            end
            %%%%% Nodal variable update
            % if W~=1
            %%%%%% Update macro velocity and density at nodes
            %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at nodes
            %     for r=1:N
            %         f_nd_eq(:,r)=eqm_h(V,U_nd(:,r),Rho_nd(1,r),T(1,r),qh,wh,Rho_r,FD);
            %     end
            %         for l=1:N;
            %             NP=NODE{l};
            %             if NP{2}==0
            %                 f_nd_eq(:,l)=node_star(NP,f_eq,0);
            %             else
            %                 f_nd_eq(:,l)=eqm_h(V,U_nd(:,l),Rho_nd(1,l),T(1,r),qh,wh,Rho_r,FD);
            %             end
            %         end
            % end
            %%%% Solving governing equation
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate collision integral over triangles
            for r=1:M
                P=CELL{r};
                FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
            end
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                %f(:,r)=f_old(:,r)+dt/P{6}*(FCOL(:,r)-FLX(:,r));
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(FCOL(:,r)-FLX(:,r));
            end
        elseif FS==1 % Collision first
            %             %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %             %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            % Collision
            for r=1:M
                f_old(:,r)=(f_old(:,r)-f_eq(:,r))*exp(-dt/Tau)+f_eq(:,r);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0; % Interior nodes
                    if W~=1 && FMP>1
                        %f_nd(:,l)=node_star(NP,f_old,0);
                    end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end
            end
            %%%%%% Update macro velocity and density at nodes
            %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %             %%%% Calculate equilibrium pdf at nodes
            %                 for r=1:N
            %                     f_nd_eq(:,r)=eqm_h(V,U_nd(:,r),Rho_nd(1,r),T(1,r),qh,wh,Rho_r,FD);
            %                 end
            % %                     for l=1:N;
            % %                         NP=NODE{l};
            % %                         if NP{2}==0
            % %                             f_nd_eq(:,l)=node_star(NP,f_eq,0);
            % %                         else
            % %                             f_nd_eq(:,l)=eqm_h(V,U_nd(:,l),Rho_nd(1,l),T(1,r),qh,wh,Rho_r,FD);
            % %                         end
            % %                     end
            %
            %             %%%% Solving governing equation
            %             %%%% Calculate flux integral for all triangles
            %             %%%% After collision pdf at cell centriod and vertices
            %             for r=1:M
            %                 f_old(:,r)=f_old(:,r)+(f_eq(:,r)-f_old(:,r))/Tau*dt;
            %             end
            %             for l=1:N;
            %                 f_nd(:,l)=f_nd(:,l)+(f_nd_eq(:,l)-f_nd(:,l))/Tau*dt;
            %             end
            
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(-FLX(:,r));
            end
        elseif FS==2 % Advection first
            %             %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %             %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    if W~=1 && FMP>1
                        %f_nd(:,l)=node_star(NP,f_old,0);
                    end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end
            end
            %%%%%% Update macro velocity and density at nodes
            %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                f_old(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(-FLX(:,r));
            end
            %             %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %             %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            % Collision
            for r=1:M
                f(:,r)=(f_old(:,r)-f_eq(:,r))*exp(-dt/Tau)+f_eq(:,r);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
        elseif FS==3 % Collision first with new method for collision
            if tt==1
                %             %%%%%% Update macro velocity and density at cell centroid
                [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
                %             %%%% Calculate equilibrium pdf at cell centroids
                for r=1:M
                    f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
                end
                %%%%f_eq_old
                f_eq_old=f_eq;
                % Collision
                for r=1:M
                    f_old(:,r)=(f_old(:,r)-f_eq(:,r))*exp(-dt/Tau)+f_eq(:,r);
                end
                %%%%%%Nodal pdf updated and boundary conditions applied
                for l=1:N;
                    NP=NODE{l};
                    if NP{2}==0; % Interior nodes
                        if W~=1 && FMP>1
                            %f_nd(:,l)=node_star(NP,f_old,0);
                        end
                    else % nodes on boundary
                        f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                    end
                end
                %%%%%% Update macro velocity and density at nodes
                %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
                %%%% Calculate flux integral for all triangles
                FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
                %%%% Calculate the current pdf at cell centroid/// See derivation in
                %%%% notes
                for r=1:M
                    P=CELL{r};
                    f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(-FLX(:,r));
                end
            else
                %             %%%%%% Update macro velocity and density at cell centroid
                [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
                %             %%%% Calculate equilibrium pdf at cell centroids
                for r=1:M
                    f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
                end
                %%%%f_eq_old
                f_eq_new=2*f_eq-f_eq_old;
                % Collision
                Grad_f_eq=(f_eq-f_eq_old)/dt;
                for r=1:M
                    f_neq(:,r)=(f_old(:,r)-f_eq(:,r)-Grad_f_eq(:,r)*Tau)*exp(-dt/Tau)+Grad_f_eq(:,r)*Tau;
                end
                f_old=f_neq+f_eq_new;
                %%%%%%Nodal pdf updated and boundary conditions applied
                for l=1:N;
                    NP=NODE{l};
                    if NP{2}==0; % Interior nodes
                        if W~=1 && FMP>1
                            %f_nd(:,l)=node_star(NP,f_old,0);
                        end
                    else % nodes on boundary
                        f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                    end
                end
                %%%%%% Update macro velocity and density at nodes
                %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
                %%%% Calculate flux integral for all triangles
                FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
                %%%% Calculate the current pdf at cell centroid/// See derivation in
                %%%% notes
                for r=1:M
                    P=CELL{r};
                    f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(-FLX(:,r));
                end
                %%%%f_eq_old
                f_eq_old=f_eq;
            end
        else
            error('The flag for splitting option is invalid or unavailable!');
        end
    elseif FT==1 % Implicit Euler, CI/AE
        if tt==1
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    % if W~=1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                    % end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end;
            end;
            %%%%% Nodal variable update
            % if W~=1
            %%%%%% Update macro velocity and density at nodes
            [Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate collision integral over triangles
            for r=1:M
                P=CELL{r};
                FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
            end
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                %f(:,r)=f_old(:,r)+dt/P{6}*(FCOL(:,r)-FLX(:,r));
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(FCOL(:,r)-FLX(:,r));
            end
            %%%%f_eq_old
            f_eq_old=f_eq;
        else
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
%                     if W~=1 && FMP>1
%                         %f_nd(:,l)=node_star(NP,f_old,0);
%                     end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end;
            end;
            %%%%% Nodal variable update
            % if W~=1
            %%%%%% Update macro velocity and density at nodes
            %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate the extrapolated f_eq_im
            f_eq_im=2*f_eq-f_eq_old;
            %%%% Calculate the current pdf at cell centroid-Implicit scheme
            for r=1:M
                P=CELL{r};
                %% 1st Algorithm, fully implicit for collision term
                f(:,r)=(f_old(:,r)-FLX(:,r)*dt/P{6}+f_eq_im(:,r)*dt/Tau)/(1+dt/Tau);
                %% 2nd Algorithm, half implicit for collision term
%                 f(:,r)=((1-dt/Tau/2)*f_old(:,r)-FLX(:,r)*dt/P{6}+(f_eq_im(:,r)+f_eq(:,r))*dt/Tau/2)/(1+dt/Tau/2);
            end
            %%%%f_eq_old
            f_eq_old=f_eq;
        end
    elseif FT==2 % Implicit Euler, CI/AI
        if tt==1
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    % if W~=1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                    % end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end;
            end;
            %%%%% Nodal variable update
            % if W~=1
            %%%%%% Update macro velocity and density at nodes
            [Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate collision integral over triangles
            for r=1:M
                P=CELL{r};
                FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
            end
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                %f(:,r)=f_old(:,r)+dt/P{6}*(FCOL(:,r)-FLX(:,r));
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(FCOL(:,r)-FLX(:,r));
            end
            %%%%f_eq_old
            f_eq_old=f_eq;
            FLX_old=FLX;
        else
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    if W~=1 && FMP>1
                        %f_nd(:,l)=node_star(NP,f_old,0);
                    end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end;
            end;
            %%%%% Nodal variable update
            % if W~=1
            %%%%%% Update macro velocity and density at nodes
            %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate the extrapolated f_eq_im and FLX_im
            f_eq_im=2*f_eq-f_eq_old;
            FLX_im=2*FLX-FLX_old;
            %%%% Calculate the current pdf at cell centroid-Implicit scheme
            for r=1:M
                P=CELL{r};
                %% 1st Algorithm, fully implicit for collision term
                f(:,r)=(f_old(:,r)-FLX_im(:,r)*dt/P{6}+f_eq_im(:,r)*dt/Tau)/(1+dt/Tau);
                %% 2nd Algorithm, half implicit for collision term
%                 f(:,r)=((1-dt/Tau/2)*f_old(:,r)-(FLX_im(:,r)+FLX(:,r))*dt/2/P{6}+(f_eq_im(:,r)+f_eq(:,r))*dt/Tau/2)/(1+dt/Tau/2);
            end
            %%%%f_eq_old
            f_eq_old=f_eq;
            FLX_old=FLX;
        end
    elseif FT==3 % AB
        if tt==1 % Fisrt time step
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N;
                NP=NODE{l};
                if NP{2}==0; % Interior nodes
                    % if W~=1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                    % end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end
            end
            %%%%%% Update macro velocity and density at nodes
            [Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Solving governing equation
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate collision integral over triangles
            for r=1:M
                P=CELL{r};
                FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
            end
            FC_old_old=FLX+FCOL;
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(FCOL(:,r)-FLX(:,r));
            end
        else % Afterward
            %%%%%% Update macro velocity and density at cell centroid
            [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
            %%%% Calculate equilibrium pdf at cell centroids
            for r=1:M
                f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
            end
            %%%%%%Nodal pdf updated and boundary conditions applied
            for l=1:N
                NP=NODE{l};
                if NP{2}==0 % Interior nodes
                    if W~=1 && FMP>1
                        %f_nd(:,l)=node_star(NP,f_old,0);
                    end
                else % nodes on boundary
                    f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
                end
            end
            %%%%%% Update macro velocity and density at nodes
            [Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
            %%%% Solving governing equation
            %%%% Calculate flux integral for all triangles
            FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
            %%%% Calculate collision integral over triangles
            for r=1:M
                P=CELL{r};
                FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
            end
            FC_old=FLX+FCOL;
            %%%% Calculate the current pdf at cell centroid/// See derivation in
            %%%% notes
            for r=1:M
                P=CELL{r};
                %% Algorithm 1
%                 f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/2/P{6}/W*(3*FC_old(:,r)-FC_old_old(:,r));
                %% Algorithm 2
                f(:,r)=f_old(:,r)+(1-W)/3/W*(f_nd_old(:,P{7})+f_nd_old(:,P{8})+f_nd_old(:,P{9})-f_nd(:,P{7})-f_nd(:,P{8})-f_nd(:,P{9}))+dt/P{6}/W*(2*FC_old(:,r)-FC_old_old(:,r));
            end
            %%%% Update
            FC_old_old=FC_old;
            f_old_old=f_old;
        end
    elseif FT==4 % RK2
        %% 1st algorithm for RK2
        %%%%%% Update macro velocity and density at cell centroid
        [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        %%%% Calculate the current pdf at cell centroid/// See derivation in
        %%%% notes
        for r=1:M
            P=CELL{r};
            f_K1(:,r)=f_old(:,r)+dt/P{6}*(FCOL(:,r)-FLX(:,r));
        end
        %%%%%% Update macro velocity and density at cell centroid
        [Rho_K1,U_K1,T_K1]=macro_h(f_K1,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq_K1(:,r)=eqm_h(V,U_K1(:,r),Rho_K1(1,r),T_K1(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%% Update macro velocity and density at nodes
        Rho_nd_K1=Rho_nd;
        U_nd_K1=U_nd;
        f_nd_K1=f_nd;
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd_K1(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_K1,f_eq_K1,f_nd_K1,U_nd_K1,U_K1,Rho_nd_K1,Rho_K1,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd_K1,U_nd_K1,T_nd_K1]=macro_h(f_nd_K1,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX_K1=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd_K1,f_K1,f_eq_K1,Rho_K1,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL_K1(:,r)=BGK(qh,P,W,Tau,f_K1(:,r),f_nd_K1(:,P{7}),f_nd_K1(:,P{8}),f_nd_K1(:,P{9}),f_eq_K1(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        %%%% Calculate the current pdf at cell centroid/// See derivation in
        %%%% notes
        for r=1:M
            P=CELL{r};
            f(:,r)=f_old(:,r)/2+f_K1(:,r)/2+dt/2/P{6}*(FCOL_K1(:,r)-FLX_K1(:,r));
        end
        
        %% 2nd algorithm for RK2
%         %%%%%% Update macro velocity and density at cell centroid
%         [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
%         %%%% Calculate equilibrium pdf at cell centroids
%         for r=1:M
%             f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
%         end
%         %%%%%%Nodal pdf updated and boundary conditions applied
%         for l=1:N
%             NP=NODE{l};
%             if NP{2}==0 % Interior nodes
%                 if W~=1 && FMP>1
%                     %f_nd(:,l)=node_star(NP,f_old,0);
%                 end
%             else % nodes on boundary
%                 f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
%             end
%         end
%         %%%%%% Update macro velocity and density at nodes
%         %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
%         %%%% Solving governing equation
%         %%%% Calculate flux integral for all triangles
%         FLX=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
%         %%%% Calculate collision integral over triangles
%         for r=1:M
%             P=CELL{r};
%             FCOL(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
%         end
%         %%%% Calculate f_K1
%         %%%% notes
%         for r=1:M
%             P=CELL{r};
%             f_K1(:,r)=f_old(:,r)+dt/P{6}/2*(FCOL(:,r)-FLX(:,r));
%         end
%         %%%%%% Update macro velocity and density at cell centroid
%         [Rho_K1,U_K1,T_K1]=macro_h(f_K1,V,Rho_r,FD);
%         %%%% Calculate equilibrium pdf at cell centroids
%         for r=1:M
%             f_eq_K1(:,r)=eqm_h(V,U_K1(:,r),Rho_K1(1,r),T_K1(1,r),qh,wh,Rho_r,FD);
%         end
%         %%%%%% Update macro velocity and density at nodes
%         Rho_nd_K1=Rho_nd;
%         U_nd_K1=U_nd;
%         f_nd_K1=f_nd;
%         %%%%%%Nodal pdf updated and boundary conditions applied
%         for l=1:N
%             NP=NODE{l};
%             if NP{2}==0 % Interior nodes
%                 if W~=1 && FMP>1
%                     %f_nd(:,l)=node_star(NP,f_old,0);
%                 end
%             else % nodes on boundary
%                 f_nd_K1(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_K1,f_eq_K1,f_nd_K1,U_nd_K1,U_K1,Rho_nd_K1,Rho_K1,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
%             end
%         end
%         %%%%%% Update macro velocity and density at nodes
%         %[Rho_nd_K1,U_nd_K1,T_nd_K1]=macro_h(f_nd_K1,V,Rho_r,FD);
%         %%%% Solving governing equation
%         %%%% Calculate flux integral for all triangles
%         FLX_K1=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd_K1,f_K1,f_eq_K1,Rho_K1,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
%         %%%% Calculate collision integral over triangles
%         for r=1:M
%             P=CELL{r};
%             FCOL_K1(:,r)=BGK(qh,P,W,Tau,f_K1(:,r),f_nd_K1(:,P{7}),f_nd_K1(:,P{8}),f_nd_K1(:,P{9}),f_eq_K1(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
%         end
%         %%%% Calculate the current pdf at cell centroid/// See derivation in
%         %%%% notes
%         for r=1:M
%             P=CELL{r};
%             f(:,r)=f_old(:,r)+dt/P{6}*(FCOL_K1(:,r)-FLX_K1(:,r));
%         end
    elseif FT==5 % RK4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Update macro velocity and density at cell centroid
        [Rho,U,T]=macro_h(f_old,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq(:,r)=eqm_h(V,U(:,r),Rho(1,r),T(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX_K1=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_old,f_eq,Rho,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL_K1(:,r)=BGK(qh,P,W,Tau,f_old(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        for r=1:M
            P=CELL{r};
            K1(:,r)=(FCOL_K1(:,r)-FLX_K1(:,r))/P{6};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_K1=f_old+dt/2*K1;
        
        [Rho_K1,U_K1,T_K1]=macro_h(f_K1,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq_K1(:,r)=eqm_h(V,U_K1(:,r),Rho_K1(1,r),T_K1(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N;
            NP=NODE{l};
            if NP{2}==0; % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_K1,f_eq_K1,f_nd,U_nd,U_K1,Rho_nd,Rho_K1,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX_K2=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_K1,f_eq_K1,Rho_K1,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL_K2(:,r)=BGK(qh,P,W,Tau,f_K1(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq_K1(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        for r=1:M
            P=CELL{r};
            K2(:,r)=(FCOL_K2(:,r)-FLX_K2(:,r))/P{6};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_K2=f_old+dt/2*K2;
        
        [Rho_K2,U_K2,T_K2]=macro_h(f_K2,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq_K2(:,r)=eqm_h(V,U_K2(:,r),Rho_K2(1,r),T_K2(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_K2,f_eq_K2,f_nd,U_nd,U_K2,Rho_nd,Rho_K2,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX_K3=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_K2,f_eq_K2,Rho_K2,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL_K3(:,r)=BGK(qh,P,W,Tau,f_K2(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq_K2(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        for r=1:M
            P=CELL{r};
            K3(:,r)=(FCOL_K3(:,r)-FLX_K3(:,r))/P{6};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_K3=f_old+dt*K3;
        
        [Rho_K3,U_K3,T_K3]=macro_h(f_K3,V,Rho_r,FD);
        %%%% Calculate equilibrium pdf at cell centroids
        for r=1:M
            f_eq_K3(:,r)=eqm_h(V,U_K3(:,r),Rho_K3(1,r),T_K3(1,r),qh,wh,Rho_r,FD);
        end
        %%%%%%Nodal pdf updated and boundary conditions applied
        for l=1:N
            NP=NODE{l};
            if NP{2}==0 % Interior nodes
                if W~=1 && FMP>1
                    %f_nd(:,l)=node_star(NP,f_old,0);
                end
            else % nodes on boundary
                f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_K3,f_eq_K3,f_nd,U_nd,U_K3,Rho_nd,Rho_K3,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
            end
        end
        %%%%%% Update macro velocity and density at nodes
        %[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
        %%%% Solving governing equation
        %%%% Calculate flux integral for all triangles
        FLX_K4=flux_pdf_h(CELL,M,NODE,N,FACE,O,f_nd,f_K3,f_eq_K3,Rho_K3,Rho_r,V,dt,wl,wlb,wh,V1,V2,FInr,FPdc,FInt,FOut,FStw,FMow,FWdp,FTVD,FPDC,FD,FM,FMP);
        %%%% Calculate collision integral over triangles
        for r=1:M
            P=CELL{r};
            FCOL_K4(:,r)=BGK(qh,P,W,Tau,f_K3(:,r),f_nd(:,P{7}),f_nd(:,P{8}),f_nd(:,P{9}),f_eq_K3(:,r),f_nd_eq(:,P{7}),f_nd_eq(:,P{8}),f_nd_eq(:,P{9}));
        end
        for r=1:M
            P=CELL{r};
            K4(:,r)=(FCOL_K4(:,r)-FLX_K4(:,r))/P{6};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FLX=(FLX_K1+FLX_K2+FLX_K3+FLX_K4)/4;
        f=f_old+dt*(K1/6+K2/3+K3/3+K4/6);
    else
        error('The flag for time-matching scheme is unavailable or invalid!');
    end
    %%%% Monitored variables
    R(tt+1)=norm(f-f_old,'fro');
    %R(tt+1)=norm(f-f_eq,'fro');
    RHO(tt+1)=0;
    for r=1:M
        RHO(tt+1)=RHO(tt+1)+Rho(1,r);
    end
    RHO(tt+1)=RHO(tt+1)/M;
    
    UR(:,tt+1)=[0;0];
    for r=1:M
        UR(:,tt+1)=UR(:,tt+1)+U(:,r);
    end
    UR(:,tt+1)=UR(:,tt+1)/M;
    
    TR(tt+1)=0;
    for r=1:M
        TR(tt+1)=TR(tt+1)+T(1,r);
    end
    TR(tt+1)=TR(tt+1)/M;
    
    FLX_c(1:qh,tt+1)=0;
    for k=1:qh
        for r=1:M
            FLX_c(k,tt+1)=FLX_c(k,tt+1)+FLX(k,r);
        end
    end
    FLX_c(1:qh,tt+1)=FLX_c(1:qh,tt+1)/M;
    
    %     FCOL_c(1:qh,tt+1)=0;
    %     for k=1:qh
    %         for r=1:M
    %             FCOL_c(k,tt+1)=FCOL_c(k,tt+1)+FCOL(k,r);
    %         end
    %     end
    %     FCOL_c(1:qh,tt+1)=FCOL_c(1:qh,tt+1)/M;
    %     if FM==0
    %         for l=1:a
    %             RHO_NDX(1+2*(l-1):2*l,tt+1)=[Rho_nd(NDX(1,l));Rho_nd(NDX(2,l))];
    %             RHO_NDY(1+2*(l-1):2*l,tt+1)=[Rho_nd(NDY(1,l));Rho_nd(NDY(2,l))];
    %             U_NDX(1+4*(l-1):4*l,tt+1)=[U_nd(:,NDX(1,l));U_nd(:,NDX(2,l))];
    %             U_NDY(1+4*(l-1):4*l,tt+1)=[U_nd(:,NDY(1,l));U_nd(:,NDY(2,l))];
    %             if FF==1 % Taylor vortex flow
    %                 U_TV_NDX(1+4*(l-1):4*l,tt+1)=[-U_0*exp(-Tau/3*dt*tt*(k1^2+k2^2))*cos(k1*X(NDX(1,l)))*sin(k2*Y(NDX(1,l)));
    %                     U_0*(k1/k2)*exp(-Tau/3*dt*tt*(k1^2+k2^2))*sin(k1*X(NDX(1,l)))*cos(k2*Y(NDX(1,l)));
    %                     -U_0*exp(-Tau/3*dt*tt*(k1^2+k2^2))*cos(k1*X(NDX(2,l)))*sin(k2*Y(NDX(2,l)));
    %                     U_0*(k1/k2)*exp(-Tau/3*dt*tt*(k1^2+k2^2))*sin(k1*X(NDX(2,l)))*cos(k2*Y(NDX(2,l)));];
    %                 U_TV_NDY(1+4*(l-1):4*l,tt+1)=[-U_0*exp(-Tau/3*dt*tt*(k1^2+k2^2))*cos(k1*X(NDY(1,l)))*sin(k2*Y(NDY(1,l)));
    %                     U_0*(k1/k2)*exp(-Tau/3*dt*tt*(k1^2+k2^2))*sin(k1*X(NDY(1,l)))*cos(k2*Y(NDY(1,l)));
    %                     -U_0*exp(-Tau/3*dt*tt*(k1^2+k2^2))*cos(k1*X(NDY(2,l)))*sin(k2*Y(NDY(2,l)));
    %                     U_0*(k1/k2)*exp(-Tau/3*dt*tt*(k1^2+k2^2))*sin(k1*X(NDY(2,l)))*cos(k2*Y(NDY(2,l)));];
    %             end
    %         end
    %         if FS==0
    %             FCOL_c(1:qh,tt+1)=0;
    %             for k=1:qh
    %                 for r=1:M
    %                     FCOL_c(k,tt+1)=FCOL_c(k,tt+1)+FCOL(k,r);
    %                 end
    %             end
    %             FCOL_c(1:qh,tt+1)=FCOL_c(1:qh,tt+1)/M;
    %
    %             for l=1:a
    %                 RHO_NDX(1+2*(l-1):2*l,tt+1)=[Rho_nd(NDX(1,l));Rho_nd(NDX(2,l))];
    %                 RHO_NDY(1+2*(l-1):2*l,tt+1)=[Rho_nd(NDY(1,l));Rho_nd(NDY(2,l))];
    %                 U_NDX(1+4*(l-1):4*l,tt+1)=[U_nd(:,NDX(1,l));U_nd(:,NDX(2,l))];
    %                 U_NDY(1+4*(l-1):4*l,tt+1)=[U_nd(:,NDY(1,l));U_nd(:,NDY(2,l))];
    %             end
    %         end
    %     end
    %%%% Update pdf at cell centroids and nodes
    f_old=f;
    f_nd_old=f_nd;
    %%%% Timer evolution
    tt=tt+1;
    TT(tt)=tt;
end
toc;
%%%%%%Nodal pdf updated and boundary conditions applied
for l=1:N
    NP=NODE{l};
    if NP{2}==0 % Interior nodes
        % if W~=1
        f_nd(:,l)=node_star(NP,f_old,0);
        % end
    else % nodes on boundary
        f_nd(:,l)=pdf_bc_h(N_L,N_H,N_I,Tau,dt,NODE,NP,f_old,f_eq,f_nd,U_nd,U,Rho_nd,Rho,V,V1,V2,Rho_r,Rho_in,Rho_out,qh,wc,wh,FD,FDC,NPdc,NInt,NOut,NStw,NMow,NWdp,NC70,NC71,NC72,NC73,NC74,NC75,NC76,NC77);
    end;
end;
[Rho_nd,U_nd,T_nd]=macro_h(f_nd,V,Rho_r,FD);
%%% Calculate the velocity megnitude in cell centroid
U_M=zeros(1,M);
for r=1:M
    P=CELL{r};
    U_M(r)=sqrt(U(1,r)^2+U(2,r)^2);
end
%%% Adding the boundary values in the matrices for display
XXX=zeros(1,M);
YYY=zeros(1,M);
%%%% coordinates
for r=1:M;
    P=CELL{r};
    Centroid=P{5};
    XXX(r)=Centroid(1,1);
    YYY(r)=Centroid(2,1);
end
for r=1:N-N_I;
    XXX(M+r)=X(r+N_I);
    YYY(M+r)=Y(r+N_I);
end
for r=1:N_I_N;
    XXX(M+N-N_I+r)=X(r);
    YYY(M+N-N_I+r)=Y(r);
end
U_plt=U;
Rho_plt=Rho;
T_plt=T;
for r=1:N-N_I;
    U_plt(:,M+r)=U_nd(:,r+N_I);
    U_M(M+r)=sqrt(U_nd(1,r+N_I)^2+U_nd(2,r+N_I)^2);
    Rho_plt(1,M+r)=Rho_nd(1,r+N_I);
    T_plt(1,M+r)=T_nd(1,r+N_I);
end
for r=1:N_I_N;
    U_plt(:,M+N-N_I+r)=U_nd(:,r);
    U_M(M+N-N_I+r)=sqrt(U_nd(1,r)^2+U_nd(2,r)^2);
    Rho_plt(1,M+N-N_I+r)=Rho_nd(1,r);
    T_plt(1,M+N-N_I+r)=T_nd(1,r);
end
Time_per_iter=toc/i