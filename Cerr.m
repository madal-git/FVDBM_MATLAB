if FTH==0
    if FF==0
        if strcmp(top,'Moving Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % Couette flow
            if FM==0 % IRT mesh
                %%% Error
                %% Hydro
                L=length(ERR(1,:));
                ErrA=abs(ERR)/U_cav(1);
                Err_sum_x=0;
                Err_sum_y=0;
                for i=1:L
                    if i==1 || i==L
                        Err_sum_x=Err_sum_x+ErrA(1,i)/2;
                        Err_sum_y=Err_sum_y+ErrA(2,i)/2;
                    else
                        Err_sum_x=Err_sum_x+ErrA(1,i);
                        Err_sum_y=Err_sum_y+ErrA(2,i);
                    end
                end
                Err_avg_x=Err_sum_x/L
                Err_avg_y=Err_sum_y/L
                
                Err_L1_x=norm(ERR(1,:),1)/norm(U_ana(1,:),1)
                Err_L2_x=norm(ERR(1,:),2)/norm(U_ana(1,:),2)
                Err_L1_y=norm(ERR(2,:),1)/norm(U_ana(1,:),1)
                Err_L2_y=norm(ERR(2,:),2)/norm(U_ana(1,:),2)
                
                %             Err_L1_x=norm(ERR(1,:)/U_cav(1),1)
                %             Err_L2_x=norm(ERR(1,:)/U_cav(1),2)
                %             Err_L1_y=norm(ERR(2,:)/U_cav(1),1)
                %             Err_L2_y=norm(ERR(2,:)/U_cav(1),2)
                
                
                Err_max_x=max(abs(ERR(1,:)))/U_cav(1)
                Err_max_y=max(abs(ERR(2,:)))/U_cav(1)
                %%% slip velocities
                for i=1:N_H
                    CXB=Y1;
                    CXT=Y2;
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXB)
                            CNXB(i)=r;
                            break;
                        end
                    end
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXT)
                            CNXT(i)=r;
                            break;
                        end
                    end
                end
                % Top wall
                for i=1:N_L
                    U_err_x_top(i)=U_cav(1)-U_nd(1,CNXT(i));
                    U_err_y_top(i)=0-U_nd(2,CNXT(i));
                end
                u_top_avg=norm(U_err_x_top,1)/N_L
                u_top_LInf=norm(U_err_x_top,Inf)
                v_top_avg=norm(U_err_y_top,1)/N_L
                v_top_LInf=norm(U_err_y_top,Inf)
                % Bottom wall
                for i=1:N_L
                    U_err_x_botm(i)=0-U_nd(1,CNXB(i));
                    U_err_y_botm(i)=0-U_nd(2,CNXB(i));
                end
                u_botm_avg=norm(U_err_x_botm,1)/N_L
                u_botm_LInf=norm(U_err_x_botm,Inf)
                v_botm_avg=norm(U_err_y_botm,1)/N_L
                v_botm_LInf=norm(U_err_y_botm,Inf)
                %% Thermal
                % Top  and bottom wall
                for i=1:N_L
                    T_err_top(i)=T_bc(1)-T_nd(1,CNXT(i));
                    T_err_botm(i)=T_bc(3)-T_nd(1,CNXB(i));
                end
                
                Err_T_top_avg=norm(T_err_top,1)/N_L;
                Err_T_botm_avg=norm(T_err_botm,1)/N_L
                %
                Err_L1_T_1=norm(T_err,1)/norm(T_ana_4err,1)
                Err_L2_T_1=norm(T_err,2)/norm(T_ana_4err,2)
                %
                %             Err_max_T=max(abs(T_err_rel))
            elseif FM==1 % Arbitrary mesh
                U_ana=U_cav/Y2*C_cut(2,1:N_cut);
                U_sim=U_cut;
                Err_L2_x=norm(U_ana(1,:)-U_sim(1,:),2)/norm(U_ana(1,:),2)
                Err_L2_y=norm(U_ana(2,:)-U_sim(2,:),2)/norm(U_ana(1,:),2)
            else
                error('Wrong flag for flow type!');
            end
        elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Periodic')~=1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')~=1% Poiseulle flow
            if FM==0 % IRT mesh
                for i=1:N_H
                    CXB=Y1;
                    CXT=Y2;
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXB)
                            CNXB(i)=r;
                            break;
                        end
                    end
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXT)
                            CNXT(i)=r;
                            break;
                        end
                    end
                end
                %             Err_u_x_out=U_err_x_out/U_ana_max_out;
                %             Err_u_y_out=U_err_y_out/U_ana_max_out;
                %             Err_rho=Rho_err/abs(Rho_l-Rho_r);
                
                Err_L1_x_out=norm(U_err_x_out,1)/norm(U_ana_x_out,1)
                Err_L2_x_out=norm(U_err_x_out,2)/norm(U_ana_x_out,2)
                Err_L1_y_out=norm(U_err_y_out,1)%/norm(U_ana_out,1)
                Err_L2_y_out=norm(U_err_y_out,2)%/norm(U_ana_out,2)
                Err_L1_rho=norm(Rho_err,1)/norm(Rho_ana,1)
                Err_L2_rho=norm(Rho_err,2)/norm(Rho_ana,2)
                
                %             Err_max_x_out=max(abs(Err_u_x_out))
                %             Err_max_y_out=max(abs(Err_u_y_out))
                %             Err_max_rho=max(abs(Err_rho))
                %%%%%%%%%%%%%%% Slip velocity
                % Inlet
                %             for i=1:N_H
                %                 ND=NODE{CNXI(i)};
                %                 BC=ND{21};
                %                 U_err_x_in_2099(i)=BC(2)-U_nd(1,CNXI(i));
                %                 U_err_y_in_2099(i)=BC(3)-U_nd(2,CNXI(i));
                %             end
                %             u_inlet_avg=norm(U_err_x_in_2099,1)/N_H
                %             u_inlet_LInf=norm(U_err_x_in_2099,Inf)
                %             v_inlet_avg=norm(U_err_y_in_2099,1)/N_H
                %             v_inlet_LInf=norm(U_err_y_in_2099,Inf)
                %             % Top wall
                %             for i=1:N_L
                %                 U_err_x_top(i)=0-U_nd(1,CNXT(i));
                %                 U_err_y_top(i)=0-U_nd(2,CNXT(i));
                %             end
                %             u_top_avg=norm(U_err_x_top,1)/N_L
                %             u_top_LInf=norm(U_err_x_top,Inf)
                %             v_top_avg=norm(U_err_y_top,1)/N_L
                %             v_top_LInf=norm(U_err_y_top,Inf)
                %             % Outlet
                %             for i=1:N_H
                %                 U_err_y_in_3199(i)=abs(0-U_nd(2,CNXO(i)));
                %             end
                %             v_outlet_avg=norm(U_err_y_in_3199,1)/N_H
                %             v_outlet_LInf=norm(U_err_y_in_3199,Inf)
                %             % Bottom wall
                %             for i=1:N_L
                %                 U_err_x_botm(i)=0-U_nd(1,CNXB(i));
                %                 U_err_y_botm(i)=0-U_nd(2,CNXB(i));
                %             end
                %             u_botm_avg=norm(U_err_x_botm,1)/N_L
                %             u_botm_LInf=norm(U_err_x_botm,Inf)
                %             v_botm_avg=norm(U_err_y_botm,1)/N_L
                %             v_botm_LInf=norm(U_err_y_botm,Inf)
                %             %%%%%%%%%%%%%%% Slip density
                %             % Inlet
                %             Rho_err_in=0;
                %             for i=1:N_H
                %                 ND=NODE{CNXI(i)};
                %                 BC=ND{21};
                %                 Rho_err_in(i)=BC(1)-Rho_nd(1,CNXI(i));
                %             end
                %             Rho_inlet_avg=norm(Rho_err_in,1)/N_H
                %             Rho_inlet_LInf=norm(Rho_err_in,Inf)
                %             % Outlet
                %             for i=1:N_H
                %                 ND=NODE{CNXO(i)};
                %                 BC=ND{21};
                %                 Rho_err_out(i)=BC(1)-Rho_nd(1,CNXO(i));
                %             end
                %             Rho_outlet_avg=norm(Rho_err_out,1)/N_H
                %             Rho_outlet_LInf=norm(Rho_err_out,Inf)
            elseif FM==1 % arbitrary mesh
                ;
            else
                error('Wrong flag for flow type!');
            end
        elseif strcmp(top,'Moving Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1% cavity flow
            %% Compare results with Ghia
            if single(U_m(1,1)*(X2-X1)/(Tau/3))==single(100) %% Re=100
                load('Ghia_Re100.mat');
            elseif single(U_m(1,1)*(X2-X1)/(Tau/3))==single(400) %% Re=400
                load('Ghia_Re400.mat');
            else
                error('Not available');
            end
            %% IRT Mesh
            if FM==0 % IRT mesh
                %% Sample x-velocity on the vertical mid plan of the numerical solution at the discrete locations of Ghia
                CX=(X2+X1)/2;
                CY=Inf;
                % Find the nodal velocity
                C_cut_y=zeros(2,N_H);
                C_cut_y(1,:)=CX;
                C_cut_y(2,:)=Y1:(Y2-Y1)/(N_H-1):Y2;
                U_cut=zeros(2,N_H);
                if mod(N_L,2)==1
                    for i=1:N_H
                        for l=1:N
                            ND=NODE{l};
                            if single(10+dis(C_cut_y(:,i),ND{3}))==single(10)
                                U_cut(:,i)=U_nd(:,l);
                                break;
                            end
                        end
                    end
                else
                    ;
                end
                % Find the velocities at Ghia's locations
                L=length(y_Ghia);
                u_Ghia_sim=zeros(L,1);
                for i=1:L
                    for j=1:N_H-1
                        if single(y_Ghia(i))>=single(C_cut_y(2,j)/(Y2-Y1)) && single(y_Ghia(i))<=single(C_cut_y(2,j+1)/(Y2-Y1))
                            break;
                        end
                    end
                    c1=abs(y_Ghia(i)-C_cut_y(2,j+1)/(Y2-Y1))/abs(C_cut_y(2,j)/(Y2-Y1)-C_cut_y(2,j+1)/(Y2-Y1));
                    c2=abs(y_Ghia(i)-C_cut_y(2,j)/(Y2-Y1))/abs(C_cut_y(2,j)/(Y2-Y1)-C_cut_y(2,j+1)/(Y2-Y1));
                    u_Ghia_sim(i)=c1*U_cut(1,j)/U_cav(1)+c2*U_cut(1,j+1)/U_cav(1);
                end
                Err_Ghia_L2_u_vs_y=norm(u_Ghia-u_Ghia_sim,2)/norm(u_Ghia,2)
                figure(1);
                hold on
                plot(u_Ghia_sim,y_Ghia);
                %% Sample x-velocity on the vertical mid plan of the numerical solution at the discrete locations of Ghia
                CX=Inf;
                CY=(Y2+Y1)/2;
                % Find the nodal velocity
                C_cut_x=zeros(2,N_L);
                C_cut_x(1,:)=X1:(X2-X1)/(N_L-1):X2;
                C_cut_x(2,:)=CY;
                V_cut=zeros(2,N_L);
                if mod(N_H,2)==1
                    for i=1:N_L
                        for l=1:N
                            ND=NODE{l};
                            if single(10+dis(C_cut_x(:,i),ND{3}))==single(10)
                                V_cut(:,i)=U_nd(:,l);
                                break;
                            end
                        end
                    end
                else
                    ;
                end
                % Find the velocities at Ghia's locations
                L=length(x_Ghia);
                v_Ghia_sim=zeros(L,1);
                for i=1:L
                    for j=1:N_L-1
                        if single(x_Ghia(i))>=single(C_cut_x(1,j)/(X2-X1)) && single(x_Ghia(i))<=single(C_cut_x(1,j+1)/(X2-X1))
                            break;
                        end
                    end
                    %                 if j==N_L-1
                    %                     if single(x_Ghia(i))==single(C_cut_x(1,j+1)/(X2-X1)) && single(x_Ghia(i))==1
                    %                         ;
                    %                     else
                    %                         error('Logic Error!');
                    %                     end
                    %                 else
                    %                     ;
                    %                 end
                    c1=abs(x_Ghia(i)-C_cut_x(1,j+1)/(X2-X1))/abs(C_cut_x(1,j)/(X2-X1)-C_cut_x(1,j+1)/(X2-X1));
                    c2=abs(x_Ghia(i)-C_cut_x(1,j)/(X2-X1))/abs(C_cut_x(1,j)/(X2-X1)-C_cut_x(1,j+1)/(X2-X1));
                    v_Ghia_sim(i)=c1*V_cut(2,j)/U_cav(1)+c2*V_cut(2,j+1)/U_cav(1);
                end
                figure(2);
                hold on
                plot(x_Ghia,v_Ghia_sim);
                Err_Ghia_L2_v_vs_x=norm(v_Ghia-v_Ghia_sim,2)/norm(v_Ghia,2)
                
                %% BVD and BDD
                %             for i=1:N_L
                %                 CXB=Y1;
                %                 CXT=Y2;
                %                 for r=1:N
                %                     if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXB)
                %                         CNXB(i)=r;
                %                         break;
                %                     end
                %                 end
                %                 for r=1:N
                %                     if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CXT)
                %                         CNXT(i)=r;
                %                         break;
                %                     end
                %                 end
                %             end
                %             for i=1:N_H
                %                 CXL=X1;
                %                 CXR=X2;
                %                 for r=1:N
                %                     if single(X(r))==single(CXL) && single(Y(r))==single(Y1+(i-1)*h)
                %                         CNXL(i)=r;
                %                         break;
                %                     end
                %                 end
                %                 for r=1:N
                %                     if single(X(r))==single(CXR) && single(Y(r))==single(Y1+(i-1)*h)
                %                         CNXR(i)=r;
                %                         break;
                %                     end
                %                 end
                %             end
                %             %%%%%%%%%%%%%%% Slip velocity
                %             % Top
                %             for i=1:N_L
                %                 ND=NODE{CNXT(i)};
                %                 BC=ND{21};
                %                 U_err_x_top(i)=BC(2)-U_nd(1,CNXT(i));
                %                 U_err_y_top(i)=BC(3)-U_nd(2,CNXT(i));
                %             end
                %             u_top_avg=norm(U_err_x_top,1)/N_L
                %             u_top_LInf=norm(U_err_x_top,Inf)
                %             v_top_avg=norm(U_err_y_top,1)/N_L
                %             v_top_LInf=norm(U_err_y_top,Inf)
                %             % Right
                %             for i=1:N_H
                %                 ND=NODE{CNXR(i)};
                %                 BC=ND{21};
                %                 U_err_x_right(i)=BC(2)-U_nd(1,CNXR(i));
                %                 U_err_y_right(i)=BC(3)-U_nd(2,CNXR(i));
                %             end
                %             u_right_avg=norm(U_err_x_right,1)/N_H
                %             u_right_LInf=norm(U_err_x_right,Inf)
                %             v_right_avg=norm(U_err_y_right,1)/N_H
                %             v_right_LInf=norm(U_err_y_right,Inf)
                %             % Bottom
                %             for i=1:N_L
                %                 ND=NODE{CNXB(i)};
                %                 BC=ND{21};
                %                 U_err_x_botm(i)=BC(2)-U_nd(1,CNXB(i));
                %                 U_err_y_botm(i)=BC(3)-U_nd(2,CNXB(i));
                %             end
                %             u_botm_avg=norm(U_err_x_botm,1)/N_L
                %             u_botm_LInf=norm(U_err_x_botm,Inf)
                %             v_botm_avg=norm(U_err_y_botm,1)/N_L
                %             v_botm_LInf=norm(U_err_y_botm,Inf)
                %             % Left
                %             for i=1:N_H
                %                 ND=NODE{CNXL(i)};
                %                 BC=ND{21};
                %                 U_err_x_left(i)=BC(2)-U_nd(1,CNXL(i));
                %                 U_err_y_left(i)=BC(3)-U_nd(2,CNXL(i));
                %             end
                %             u_left_avg=norm(U_err_x_left,1)/N_H
                %             u_left_LInf=norm(U_err_x_left,Inf)
                %             v_left_avg=norm(U_err_y_left,1)/N_H
                %             v_left_LInf=norm(U_err_y_left,Inf)
            elseif FM==1 % Arbitrary mesh
                %             figure(1);
                %             plot(u_Ghia,y_Ghia,U_cut(1,1:N_H)/U_cav(1),C_cut_y(2,1:N_cut_x)/(Y2-Y1));
                %             figure(2);
                %             plot(x_Ghia,v_Ghia,C_cut_x(1,1:N_cut_y)/(X2-X1),V_cut(2,1:N_H)/U_cav(1));
                %% Sample x-velocity on the vertical mid plan of the numerical solution at the discrete locations of Ghia
                L=length(y_Ghia);
                u_Ghia_sim=zeros(L,1);
                for i=1:L
                    for j=1:N_H-1
                        if single(y_Ghia(i))>=single(C_cut_y(2,j)/(Y2-Y1)) && single(y_Ghia(i))<=single(C_cut_y(2,j+1)/(Y2-Y1))
                            break;
                        end
                    end
                    %                 if j==N_H-1
                    %                     if single(y_Ghia(i))==single(C_cut_y(2,j+1)/(Y2-Y1)) && single(y_Ghia(i))==1
                    %                         ;
                    %                     else
                    %                         error('Logic Error!');
                    %                     end
                    %                 else
                    %                     ;
                    %                 end
                    c1=abs(y_Ghia(i)-C_cut_y(2,j+1)/(Y2-Y1))/abs(C_cut_y(2,j)/(Y2-Y1)-C_cut_y(2,j+1)/(Y2-Y1));
                    c2=abs(y_Ghia(i)-C_cut_y(2,j)/(Y2-Y1))/abs(C_cut_y(2,j)/(Y2-Y1)-C_cut_y(2,j+1)/(Y2-Y1));
                    u_Ghia_sim(i)=c1*U_cut(1,j)/U_cav(1)+c2*U_cut(1,j+1)/U_cav(1);
                end
                %             figure(1);
                %             hold on
                %             plot(u_Ghia_sim,y_Ghia);
                Err_Ghia_L2_u_vs_y=norm(u_Ghia-u_Ghia_sim,2)/norm(u_Ghia,2)
                %% Sample x-velocity on the vertical mid plan of the numerical solution at the discrete locations of Ghia
                L=length(x_Ghia);
                v_Ghia_sim=zeros(L,1);
                for i=1:L
                    for j=1:N_L-1
                        if single(x_Ghia(i))>=single(C_cut_x(1,j)/(X2-X1)) && single(x_Ghia(i))<=single(C_cut_x(1,j+1)/(X2-X1))
                            break;
                        end
                    end
                    %                 if j==N_L-1
                    %                     if single(x_Ghia(i))==single(C_cut_x(1,j+1)/(X2-X1)) && single(x_Ghia(i))==1
                    %                         ;
                    %                     else
                    %                         error('Logic Error!');
                    %                     end
                    %                 else
                    %                     ;
                    %                 end
                    c1=abs(x_Ghia(i)-C_cut_x(1,j+1)/(X2-X1))/abs(C_cut_x(1,j)/(X2-X1)-C_cut_x(1,j+1)/(X2-X1));
                    c2=abs(x_Ghia(i)-C_cut_x(1,j)/(X2-X1))/abs(C_cut_x(1,j)/(X2-X1)-C_cut_x(1,j+1)/(X2-X1));
                    v_Ghia_sim(i)=c1*V_cut(2,j)/U_cav(1)+c2*V_cut(2,j+1)/U_cav(1);
                end
                %             figure(2);
                %             hold on
                %             plot(x_Ghia,v_Ghia_sim);
                Err_Ghia_L2_v_vs_x=norm(v_Ghia-v_Ghia_sim,2)/norm(v_Ghia,2)
            else
                error('Wrong flag for mesh type!');
            end
        elseif strcmp(left,'Pressure Inlet')==1 && strcmp(right,'Pressure Outlet')==1  % Pressure driven flow
            Err_L1_rho=norm(Rho_err,1)/norm(Rho_ana,1)
            Err_L2_rho=norm(Rho_err,2)/norm(Rho_ana,2)
            %% BDD
            for i=1:N_H
                CXL=X1;
                CXR=X2;
                for r=1:N
                    if single(X(r))==single(CXL) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXI(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXR) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXO(i)=r;
                        break;
                    end
                end
            end
            % Inlet
            Rho_err_in=0;
            for i=1:N_H
                ND=NODE{CNXI(i)};
                BC=ND{21};
                Rho_err_in(i)=BC(1)-Rho_nd(1,CNXI(i));
            end
            Rho_inlet_avg=norm(Rho_err_in,1)/N_H
            Rho_inlet_LInf=norm(Rho_err_in,Inf)
            % Outlet
            for i=1:N_H
                ND=NODE{CNXO(i)};
                BC=ND{21};
                Rho_err_out(i)=BC(1)-Rho_nd(1,CNXO(i));
            end
            Rho_outlet_avg=norm(Rho_err_out,1)/N_H
            Rho_outlet_LInf=norm(Rho_err_out,Inf)
        elseif strcmp(right,'Periodic')~=1 && strcmp(left,'Periodic')~=1% Other channel flow
            for i=1:N_L
                CXB=Y1;
                CXT=Y2;
                for r=1:N
                    if single(X(r))==single(X1+(i-1)*dx) && single(Y(r))==single(CXB)
                        CNXB(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(X1+(i-1)*dx) && single(Y(r))==single(CXT)
                        CNXT(i)=r;
                        break;
                    end
                end
            end
            for i=1:N_H
                CXL=X1;
                CXR=X2;
                for r=1:N
                    if single(X(r))==single(CXL) && single(Y(r))==single(Y1+(i-1)*dy)
                        CNXL(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXR) && single(Y(r))==single(Y1+(i-1)*dy)
                        CNXR(i)=r;
                        break;
                    end
                end
            end
            %% BVD
            % inlet/left
            for i=1:N_H
                ND=NODE{CNXL(i)};
                BC=ND{21};
                U_err_x_left(i)=BC(2)-U_nd(1,CNXL(i));
                U_err_y_left(i)=BC(3)-U_nd(2,CNXL(i));
            end
            u_left_avg=norm(U_err_x_left,1)/N_H
            u_left_LInf=norm(U_err_x_left,Inf)
            v_left_avg=norm(U_err_y_left,1)/N_H
            v_left_LInf=norm(U_err_y_left,Inf)
            
            % Immersed boundary
            if N_I_N~=0 % There is immersed boundary
                angle=zeros(1,N_I_N+1);
                angle_temp=zeros(1,N_I_N+1);
                U_err_x_immersed=zeros(1,N_I_N+1);
                U_err_y_immersed=zeros(1,N_I_N+1);
                U_err_x_immersed_temp=zeros(1,N_I_N+1);
                U_err_y_immersed_temp=zeros(1,N_I_N+1);
                for r=1:N_I_N
                    ND=NODE{r};
                    if ND{2}>0
                        error('The node is not on the immersed boundary!');
                    end
                    BC=ND{21};
                    U_err_x_immersed(r)=BC(2)-U_nd(1,r);
                    U_err_y_immersed(r)=BC(3)-U_nd(2,r);
                    % angle
                    dis_horizontal=(X(r)-3)/0.3;
                    dis_vertical=(Y(r)-2.7)/0.3;
                    angle1=acos(dis_horizontal)/pi*180;
                    angle2=asin(dis_vertical)/pi*180;
                    if angle2>=0 % from 0 to 180
                        if angle1<=90 && angle2<=90 % from 0 to 90
                            angle(r)=angle2;
                        else % from 90 to 180
                            angle(r)=angle1;
                        end
                    else % from 180 to 360
                        if angle1>=90 % from 180 to 270
                            angle(r)=180+abs(angle2);
                        else % from 270 to 360
                            angle(r)=360-abs(angle2);
                        end
                    end
                end
            end
            front_stag_found=0;
            for s=1:r
                if single(angle(s))==single(180)
                    front_stag_found=1;
                    break;
                end
            end
            if front_stag_found
                % Moving the front stagnation point to the first position of the
                % data and rearrange all data input
                angle_temp(1:r-s+1)=angle(s:r)-180;
                angle_temp(r-s+2:r)=angle(1:s-1)+180;
                U_err_x_immersed_temp(1:r-s+1)=U_err_x_immersed(s:r);
                U_err_x_immersed_temp(r-s+2:r)=U_err_x_immersed(1:s-1);
                U_err_y_immersed_temp(1:r-s+1)=U_err_y_immersed(s:r);
                U_err_y_immersed_temp(r-s+2:r)=U_err_y_immersed(1:s-1);
                angle=angle_temp;
                U_err_x_immersed=U_err_x_immersed_temp;
                U_err_y_immersed=U_err_y_immersed_temp;
            end
            % filling the last position
            angle(r+1)=360+angle(1);
            U_err_x_immersed(r+1)=U_err_x_immersed(1);
            U_err_y_immersed(r+1)=U_err_y_immersed(1);
            figure
            plot(angle,U_err_x_immersed,angle,U_err_y_immersed);
            u_immersed_avg=norm(U_err_x_immersed,1)/(N_I_N+1)
            u_immersed_LInf=norm(U_err_x_immersed,Inf)
            v_immersed_avg=norm(U_err_y_immersed,1)/(N_I_N+1)
            v_immersed_LInf=norm(U_err_y_immersed,Inf)
            %% Force evaluation
            [F_p_total,F_m_n_total,F_m_s_total]=force_boundary_rec('Immersed',X1,X2,Y1,Y2,NODE,CELL,FACE,Rho,U,Rho_nd,U_nd,V,Tau,2)
            F_total=F_p_total+F_m_n_total+F_m_s_total
            %         Rho_avg=0;
            %         edge_found_counter=0;
            %         for r=1:M
            %             P=CELL{r};
            %             boundary_edge_found=0;
            %             for i=1:3
            %                 if single(P{19+(i-1)})<0
            %                     boundary_edge_found=1;
            %                     break;
            %                 end
            %             end
            %             if boundary_edge_found==1
            %                 edge_found_counter=edge_found_counter+1;
            %                 tri(edge_found_counter)=r;
            %                 edge(edge_found_counter)=i;
            %             end
            %         end
            %         for s=1:edge_found_counter
            %             Rho_avg=Rho_avg+Rho(tri(s));
            %         end
            %         Rho_avg=Rho_avg/edge_found_counter;
            Rho_avg=norm(Rho,1)/M;
            C_D=F_total(1)/Rho_avg/0.3/(U_in(1,4))^2
            C_L=F_total(2)/Rho_avg/0.3/(U_in(1,4))^2
            
            %% Pressure coefficient
            Rho_inf=zeros(1,N_H);
            for i=1:N_H
                ND=NODE{CNXR(i)};
                Rho_inf(i)=Rho_nd(1,CNXR(i));
            end
            Rho_inf=norm(Rho_inf,1)/N_H;
            if q==7
                p_inf=Rho_inf/4;
            elseif q==9
                p_inf=Rho_inf/3;
            elseif q==13
                p_inf=Rho_inf/2;
            else
                error('Other lattice is not available!');
            end
            C_P=(p_cir-p_inf)/(0.5*Rho_inf*(U_in(1,4))^2);
            figure
            plot(Angle,C_P);
        else
            msg=['The current type of flow cannot be plotted!'];
            disp(msg);
        end
    elseif FF==1 % TGV Flow
        Err_L1_vert_x=norm(ERR_vert(1,:),1)/norm(U_ana_vert(1,:),1)
        Err_L2_vert_x=norm(ERR_vert(1,:),2)/norm(U_ana_vert(1,:),2)
        Err_L1_vert_y=norm(ERR_vert(2,:),1)/U_0
        Err_L2_vert_y=norm(ERR_vert(2,:),2)/U_0
        
        Err_L1_hori_x=norm(ERR_hori(1,:),1)/U_0
        Err_L2_hori_x=norm(ERR_hori(1,:),2)/U_0
        Err_L1_hori_y=norm(ERR_hori(2,:),1)/norm(U_ana_hori(2,:),1)
        Err_L2_hori_y=norm(ERR_hori(2,:),2)/norm(U_ana_hori(2,:),2)
        
    else
        error('Other flow types are not available!');
    end
elseif FTH==1
    if FF==0
        if strcmp(top,'Periodic')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Periodic')==1 && strcmp(left,'Periodic')==1 % All periodic
            if N_H==19 && N_L==19
                load('T_ana_sampled_18by18.mat');
            elseif N_H==37 && N_L==37
                load('T_ana_sampled_36by36.mat');
            elseif N_H==73 && N_L==73
                load('T_ana_sampled_72by72.mat');
            else
                error('Check the nodal spacings on the boundaries');
            end
            T_simu=zeros(N_H,N_L);
            if FM==0 %IRT mesh
                if mod(N_L,2)==1
                    for j=1:N_L
                        C_cut_y=zeros(2,N_H);
                        C_cut_y(1,:)=X1+(X2-X1)/(N_L-1)*(j-1);
                        C_cut_y(2,:)=Y1:(Y2-Y1)/(N_H-1):Y2;
                        for i=1:N_H
                            for l=1:N
                                ND=NODE{l};
                                if single(10+dis(C_cut_y(:,i),ND{3}))==single(10)
                                    T_simu(j,i)=T_nd(1,l);
                                    break;
                                end
                            end
                        end
                    end
                else
                    error('Temporarily not available!');
                end
            elseif FM==1 % General Mesh
                for i=1:N_H
                    if i==1
                        T_simu(:,N_H)=[T_nd(1,N),T_nd(1,N_I+1:N_I+N_L-1)];
                    elseif i==N_H
                        T_simu(:,1)=fliplr(T_nd(1,N-(N_H-1)-(N_L-1):N-(N_H-1)));
                    else
                        cut=[Inf;Y2-dy*(i-1)];
                        [C_cut_x,T_simu(:,N_H-(i-1))]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,T,T_nd,CELL,M,NODE,N);
                    end
                end
            else
                error('Wrong flag for mesh type!');
            end
            surf(T_simu);
            hold on
            surf(T_ana_sampled)
            T_err_L2=norm(T_simu-T_ana_sampled,2)/norm(T_ana_sampled,2)
        elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1% 1D conduction
            %% Analytical solution
            T_ana=zeros(1,N_H);
            y_ana=0:(Y2-Y1)/(N_H-1):(Y2-Y1);
            if length(y_ana)~=N_H
                error('Error in vector length!');
            end
            for i=1:N_H
                T_ana(i)=sin(pi*y_ana(1,i)/(Y2-Y1))*exp(-Tau_t/3*2*(pi/(Y2-Y1))^2*Time_stop)+T_ref;
            end
            plot(y_ana,T_ana)
            if FM==0 % IRT mesh
                %% Numerical solution
                CX=(X2+X1)/2;
                CY=Inf;
                % Find the nodal temperature
                C_cut_y=zeros(2,N_H);
                C_cut_y(1,:)=CX;
                C_cut_y(2,:)=Y1:(Y2-Y1)/(N_H-1):Y2;
                T_cut=zeros(1,N_H);
                if mod(N_L,2)==1
                    for i=1:N_H
                        for l=1:N
                            ND=NODE{l};
                            if single(10+dis(C_cut_y(:,i),ND{3}))==single(10)
                                T_cut(1,i)=T_nd(1,l);
                                break;
                            end
                        end
                    end
                else
                    error('Temporarily not available!');
                end
                hold on
                plot(C_cut_y(2,:),T_cut(1,:))
                %% Calculate the error
                T_Err_L2=norm(T_cut-T_ana,2)/norm(T_ana,2)
            elseif FM==1 % General Mesh
                %% Numerical solution
                CX=(X2+X1)/2;
                CY=Inf;
                cut=[CX;CY];
                % Find the nodal temperature
                [C_cut_y,T_cut]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,T,T_nd,CELL,M,NODE,N);
                hold on
                plot(C_cut_y(2,:),T_cut)
                %% Calculate the error
                T_Err_L2=norm(T_cut-T_ana,2)/norm(T_ana,2)
            else
                error('Wrong flag for mesh type!');
            end
        elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1% 2D conduction with all periodic bc
            if N_H==19 && N_L==19
                load('T_ana_18by18.mat');
            elseif N_H==37 && N_L==37
                load('T_ana_36by36.mat');
            elseif N_H==73 && N_L==73
                load('T_ana_72by72.mat');
            else
                error('Check the nodal spacings on the boundaries');
            end
            T_simu=zeros(N_H,N_L);
            if FM==0 %IRT mesh
                if mod(N_L,2)==1
                    for j=1:N_L
                        C_cut_y=zeros(2,N_H);
                        C_cut_y(1,:)=X1+(X2-X1)/(N_L-1)*(j-1);
                        C_cut_y(2,:)=Y1:(Y2-Y1)/(N_H-1):Y2;
                        for i=1:N_H
                            for l=1:N
                                ND=NODE{l};
                                if single(10+dis(C_cut_y(:,i),ND{3}))==single(10)
                                    T_simu(j,i)=T_nd(1,l);
                                    break;
                                end
                            end
                        end
                    end
                else
                    error('Temporarily not available!');
                end
            elseif FM==1 % General Mesh
                for i=1:N_H
                    if i==1
                        T_simu(:,N_H)=[T_nd(1,N),T_nd(1,N_I+1:N_I+N_L-1)];
                    elseif i==N_H
                        T_simu(:,1)=fliplr(T_nd(1,N-(N_H-1)-(N_L-1):N-(N_H-1)));
                    else
                        cut=[Inf;Y2-dy*(i-1)];
                        [C_cut_x,T_simu(:,N_H-(i-1))]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,T,T_nd,CELL,M,NODE,N);
                    end
                end
            else
                error('Wrong flag for mesh type!');
            end
            surf(T_simu);
            hold on
            surf(T_ana_sampled)
            T_err_L2=norm(T_simu-T_ana_sampled,2)/norm(T_ana_sampled,2)
        elseif strcmp(top,'Fully Developed')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Fully Developed')==1 && strcmp(left,'Stationary Wall')==1 % 1D conduction
            T_ana=zeros(1,N_L);
            x_ana=0:(X2-X1)/(N_L-1):(X2-X1);
            if length(x_ana)~=N_L
                error('Error in vector length!');
            end
            for i=1:N_L
                T_ana(i)=T_bc(4)+(T_bc(4)-T_bc(2))/(X1-X2)*(x_ana(i)-X1);
            end
            plot(x_ana,T_ana)
            if FM==0 % IRT mesh
                %% Numerical solution
                CX=Inf;
                CY=(Y2+Y1)/2;
                % Find the nodal temperature
                C_cut_x=zeros(2,N_L);
                C_cut_x(1,:)=X1:(X2-X1)/(N_L-1):X2;
                C_cut_x(2,:)=CY;
                T_cut=zeros(1,N_L);
                if mod(N_H,2)==1
                    for i=1:N_L
                        for l=1:N
                            ND=NODE{l};
                            if single(10+dis(C_cut_x(:,i),ND{3}))==single(10)
                                T_cut(1,i)=T_nd(1,l);
                                break;
                            end
                        end
                    end
                else
                    error('Temporarily not available!');
                end
                hold on
                plot(C_cut_x(2,:),T_cut(1,:))
                %% Calculate the error
                T_Err_L2=norm(T_cut-T_ana,2)/norm(T_ana,2)
            elseif FM==1 % General Mesh
                %% Numerical solution
                CX=Inf;
                CY=(Y2+Y1)/2;
                cut=[CX;CY];
                % Find the nodal temperature
                [C_cut_x,T_cut]=slice(cut,N_L,N_H,N_I,X1,X2,Y1,Y2,T,T_nd,CELL,M,NODE,N);
                hold on
                plot(C_cut_x(1,:),T_cut)
                %% Calculate the error
                T_Err_L2=norm(T_cut-T_ana,2)/norm(T_ana,2)
            else
                error('Wrong flag for mesh type!');
            end
        else
            error('Other types of conduction problems are not available!');
        end
    elseif FF==1
        ;
    else
        error('Other flow types are not available!');
    end
else
    error('Wrong flag for thermal');
end