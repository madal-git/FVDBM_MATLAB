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
            L_T=length(T_err);
            
            Err_avg_T=abs(sum(T_err)/L_T)
            
            Err_L1_T_1=norm(T_err,1)/norm(T_ana_4err,1)
            Err_L2_T_1=norm(T_err,2)/norm(T_ana_4err,2)
            
            Err_max_T=max(abs(T_err_rel))
        elseif FM==1 % Arbitrary mesh
            ;
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
            Err_u_x_out=U_err_x_out/U_ana_max_out;
            Err_u_y_out=U_err_y_out/U_ana_max_out;
            Err_rho=Rho_err/abs(Rho_l-Rho_r);
            
            Err_L1_x_out=norm(U_err_x_out,1)/norm(U_ana_out,1)
            Err_L2_x_out=norm(U_err_x_out,2)/norm(U_ana_out,2)
            Err_L1_y_out=norm(U_err_y_out,1)%/norm(U_ana_out,1)
            Err_L2_y_out=norm(U_err_y_out,2)%/norm(U_ana_out,2)
            Err_L1_rho=norm(Rho_err,1)/norm(Rho_ana,1)
            Err_L2_rho=norm(Rho_err,2)/norm(Rho_ana,2)
            
            Err_max_x_out=max(abs(Err_u_x_out))
            Err_max_y_out=max(abs(Err_u_y_out))
            Err_max_rho=max(abs(Err_rho))
            %%%%%%%%%%%%%%% Slip velocity
            % Inlet
            for i=1:N_H
                ND=NODE{CNXI(i)};
                BC=ND{21};
                U_err_x_in_2099(i)=BC(2)-U_nd(1,CNXI(i));
                U_err_y_in_2099(i)=BC(3)-U_nd(2,CNXI(i));
            end
            u_inlet_avg=norm(U_err_x_in_2099,1)/N_H
            u_inlet_LInf=norm(U_err_x_in_2099,Inf)
            v_inlet_avg=norm(U_err_y_in_2099,1)/N_H
            v_inlet_LInf=norm(U_err_y_in_2099,Inf)
            % Top wall
            for i=1:N_L
                U_err_x_top(i)=0-U_nd(1,CNXT(i));
                U_err_y_top(i)=0-U_nd(2,CNXT(i));
            end
            u_top_avg=norm(U_err_x_top,1)/N_L
            u_top_LInf=norm(U_err_x_top,Inf)
            v_top_avg=norm(U_err_y_top,1)/N_L
            v_top_LInf=norm(U_err_y_top,Inf)
            % Outlet
            for i=1:N_H
                U_err_y_in_3199(i)=abs(0-U_nd(2,CNXO(i)));
            end
            v_outlet_avg=norm(U_err_y_in_3199,1)/N_H
            v_outlet_LInf=norm(U_err_y_in_3199,Inf)
            % Bottom wall
            for i=1:N_L
                U_err_x_botm(i)=0-U_nd(1,CNXB(i));
                U_err_y_botm(i)=0-U_nd(2,CNXB(i));
            end
            u_botm_avg=norm(U_err_x_botm,1)/N_L
            u_botm_LInf=norm(U_err_x_botm,Inf)
            v_botm_avg=norm(U_err_y_botm,1)/N_L
            v_botm_LInf=norm(U_err_y_botm,Inf)
            %%%%%%%%%%%%%%% Slip density
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
        elseif FM==1 % arbitrary mesh
            ;
        else
            error('Wrong flag for flow type!');
        end
    elseif strcmp(top,'Moving Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1% cavity flow
        if FM==0 % IRT mesh
            for i=1:N_L
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
            for i=1:N_H
                CXL=X1;
                CXR=X2;
                for r=1:N
                    if single(X(r))==single(CXL) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXL(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXR) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXR(i)=r;
                        break;
                    end
                end
            end
            %%%%%%%%%%%%%%% Slip velocity
            % Top
            for i=1:N_L
                ND=NODE{CNXT(i)};
                BC=ND{21};
                U_err_x_top(i)=BC(2)-U_nd(1,CNXT(i));
                U_err_y_top(i)=BC(3)-U_nd(2,CNXT(i));
            end
            u_top_avg=norm(U_err_x_top,1)/N_L
            u_top_LInf=norm(U_err_x_top,Inf)
            v_top_avg=norm(U_err_y_top,1)/N_L
            v_top_LInf=norm(U_err_y_top,Inf)
            % Right
            for i=1:N_H
                ND=NODE{CNXR(i)};
                BC=ND{21};
                U_err_x_right(i)=BC(2)-U_nd(1,CNXR(i));
                U_err_y_right(i)=BC(3)-U_nd(2,CNXR(i));
            end
            u_right_avg=norm(U_err_x_right,1)/N_H
            u_right_LInf=norm(U_err_x_right,Inf)
            v_right_avg=norm(U_err_y_right,1)/N_H
            v_right_LInf=norm(U_err_y_right,Inf)
            % Bottom
            for i=1:N_L
                ND=NODE{CNXB(i)};
                BC=ND{21};
                U_err_x_botm(i)=BC(2)-U_nd(1,CNXB(i));
                U_err_y_botm(i)=BC(3)-U_nd(2,CNXB(i));
            end
            u_botm_avg=norm(U_err_x_botm,1)/N_L
            u_botm_LInf=norm(U_err_x_botm,Inf)
            v_botm_avg=norm(U_err_y_botm,1)/N_L
            v_botm_LInf=norm(U_err_y_botm,Inf)
            % Left
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
        elseif FM==1 % Arbitrary mesh
            ;
        else
            error('Wrong flag for mesh type!');
        end
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
        [F_p_total,F_m_n_total,F_m_s_total]=force_boundary_rec('Immersed',X1,X2,Y1,Y2,NODE,CELL,UWD,Rho,U,Rho_nd,U_nd,V,Tau,2)
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
    else
        msg=['The current type of flow cannot be plotted!'];
        disp(msg);
    end
else
    error('Other flow types are not available!');
end