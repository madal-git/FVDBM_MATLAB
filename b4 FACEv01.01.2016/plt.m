figure(4);
xlin=linspace(X1,X2,100);
ylin=linspace(Y1,Y2,100);
[Xx,Yy]=meshgrid(xlin,ylin);
Z=griddata(XXX,YYY,Rho_plt,Xx,Yy);
surf(Xx,Yy,Z);
axis equal tight;
hold on;
plot3(XXX,YYY,Rho_plt,'.','MarkerSize',15);
figure(5);
contourf(Xx,Yy,Z,100);
axis equal tight;
figure(6);
Z=griddata(XXX,YYY,U_M,Xx,Yy);
contourf(Xx,Yy,Z,100);
axis equal tight;
Z=griddata(X,Y,Rho_nd,Xx,Yy);
figure(7);
contourf(Xx,Yy,Z,300);
axis equal tight;
for r=1:N
    U_nd_M(r)=sqrt(U_nd(1,r)^2+U_nd(2,r)^2);
end
figure(8);
Z=griddata(X,Y,U_nd_M,Xx,Yy);
contourf(Xx,Yy,Z,300);
axis equal tight;
if qh==37
    figure(9)
    Z=griddata(XXX,YYY,T_plt,Xx,Yy);
    contourf(Xx,Yy,Z,300);
    axis equal tight;
end
for r=1:M
    P=CELL{r};
    Centroid=P{5};
    XX(r)=Centroid(1,1);
    YY(r)=Centroid(2,1);
end
figure(10);
quiver(XX,YY,U_plt(1,1:M),U(2,1:M),10);
axis equal tight  
figure(11);
plot(TT,log10(R));
figure(12);
plot(TT,RHO);
figure(13);
plot(TT,UR(1,:),TT,UR(2,:));
figure(14);
plot(TT,TR);
figure(15);
for k=1:qh
    plot(TT,FLX_c(k,:));
    hold on
end
hold off

% if FS==0
% figure(14);
%     for k=1:qh
%         plot(TT,FCOL_c(k,:));
%         hold on
%     end
% hold off
% end

% figure(15);
% for k=1:2*a
%     plot(TT,RHO_NDX(k,:));
%     hold on
% end
% hold off
% 
% figure(16);
% for k=1:a
%     plot(TT,(RHO_NDX(1+2*(k-1),:)-RHO_NDX(2*k,:)));
%     hold on
% end
% hold off
% 
% figure(17);
% for k=1:2*a
%     plot(TT,RHO_NDY(k,:));
%     hold on
% end
% hold off
% 
% figure(18);
% for k=1:a
%     plot(TT,(RHO_NDY(1+2*(k-1),:)-RHO_NDY(2*k,:)));
%     hold on
% end
% hold off
% 
% if FF==1% Taylor vortex flow
%     figure(19);
%     for k=1:2*a  % 2*a for showing one point; 4*a for showing 2 points
%         plot(TT,log10(U_NDX(k,:)));
%         hold on
%         plot(TT,log10(U_TV_NDX(k,:)));
%         hold on
%     end
%     hold off
% end
% 
% figure(20);
% for k=1:a
%     plot(TT,(abs(U_NDX(1+4*(k-1),:))-abs(U_NDX(3+4*(k-1),:))),TT,(abs(U_NDX(2+4*(k-1),:))-abs(U_NDX(4+4*(k-1),:))));
%     hold on
% end
% hold off

if FF==0 %Normal flow
    if strcmp(top,'Moving Wall')==1 && strcmp(right,'Periodic')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')==1 % Couette flow
        if FM==0 % IRT mesh
            %% Velocity profile
            U_cav=U_m(:,1);
            for i=1:N_H
                CX=floor(N_L/2)*h;
                CXL=X1;
                %CXL=floor(N_L/2)*h-h;
                CXR=X2;
                %CXR=floor(N_L/2)*h+h;
                for r=1:N
                    if single(X(r))==single(CX) && single(Y(r))==single((i-1)*h)
                        CNX(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXL) && single(Y(r))==single((i-1)*h)
                        CNXL(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXR) && single(Y(r))==single((i-1)*h)
                        CNXR(i)=r;
                        break;
                    end
                end
            end
            % X velocity
            figure;
            plot(1/Y2*Y(CNX(1:N_H)),Y(CNX(1:N_H))/(Y2-Y1),U_nd(1,CNX(1:N_H))/U_cav(1),Y(CNX(1:N_H))/(Y2-Y1));
            title('Velocity profile of Couette flow');
            xlabel('(u-u_t)/u_t');
            ylabel('y/H');
            legend('Analytic','Simulation');
            % Y velocity
            figure;
            AnaCout=zeros(1,N_H);
            plot(Y(CNX(1:N_H)),U_nd(2,CNX(1:N_H)),Y(CNX(1:N_H)),AnaCout);
            %%%%%%%% Error
            ERR=0;
            for i=1:N_H
                U_ana(1,i)=U_cav(1)/(Y2-Y1)*Y(CNX(i));
                U_ana(2,i)=0;
                %ERR(1,i)=U_ana(1,i)-(U_nd(1,CNXL(i))+U_nd(1,CNXR(i)))/2;
                ERR(1,i)=U_ana(1,i)-U_nd(1,CNX(i));
                ERR(2,i)=U_ana(2,i)-U_nd(2,CNX(i));
                %ERR(2,i)=U_nd(2,CNX(i))-0;
            end
            figure
            plot(Y(CNX(1:N_H)),ERR(1,1:N_H)/U_cav(1));
            %% Temperature profile
            if qh==37
                %% Analytic solution
                Step_Y=1/1000;
                Y_ana=0:Step_Y:1;
                Pr=1; % For BGK and single set of PDF for both hydrodynamics and thermodynamics
                Cp=2; % Cp=(D/2)+1
                Ec=U_m(1,1)^2/Cp/(T_bc(1)-T_bc(3));
                T_ana=Y_ana+Pr*Ec/2*(Y_ana.*(1-Y_ana));
                %% plot
                figure;
                plot(T_ana,Y_ana,(T_nd(1,CNX(1:N_H))-T_bc(3))/(T_bc(1)-T_bc(3)),Y(CNX(1:N_H))/(Y2-Y1));
                title(['Temperature profile of thermal Couette flow at Br=Pr*Ec=', num2str(Pr*Ec)]);
%                 xlabel('$$\theta-\theta_b\over\theta_t-\theta_b$$','Interpreter', 'Latex');
                xlabel('(\theta-\theta_b)/(\theta_t-\theta_b)');
                ylabel('y/H');
                legend('Analytic','Simulation');
                %% calculate error and plot
                Y_ana_4err=zeros(1,N_H);
                for i=1:N_H
                    Y_ana_4err(i)=Y(CNX(i))/(Y2-Y1);
                end
                T_ana_4err=Y_ana_4err+Pr*Ec/2*(Y_ana_4err.*(1-Y_ana_4err));
                T_err=zeros(1,N_H);
                for i=1:N_H
                    T_err(i)=(T_nd(1,CNX(i))-T_bc(3))/(T_bc(1)-T_bc(3))-T_ana_4err(i);
                end
                T_err_rel=T_err./T_ana_4err;
                T_err_rel(1)=0; % since T_ana_4err(1)=0, T_err_rel(1)=inf. So, force it to be zero
                figure
                plot(T_err_rel,Y(CNX(1:N_H))/(Y2-Y1));
                title('Relative error of temperature');
                xlabel('(\theta_s_i_m-\theta_a_n_a)/\theta_a_n_a');
                ylabel('y/H');
            else
                disp('Temperature profile for thermal Couette flow is temporarilly not available!');
            end
        elseif FM==1 % Arbitrary mesh
            U_cav=U_m(:,1);
            CX=(X2+X1)/2;
            CY=Inf;
            cut=[CX;CY];
            [C_cut,U_cut]=slice(cut,N_L,N_H,X1,X2,Y1,Y2,U,U_nd,CELL,M);
            N_cut=length(U_cut);
            plot(C_cut(2,1:N_cut),U_cut(1,1:N_H),C_cut(2,1:N_cut),U_cav(1)/Y2*C_cut(2,1:N_cut));
        else
            ;
        end
    elseif strcmp(top,'Stationary Wall')==1 && strcmp(right,'Periodic')~=1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Periodic')~=1% Poiseulle flow
        if FM==0 % IRT mesh
            for i=1:N_H
                CXI=X1;
                CXO=X2;
                for r=1:N
                    if single(X(r))==single(CXI) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXI(i)=r;
                        break;
                    end
                end
                for r=1:N
                    if single(X(r))==single(CXO) && single(Y(r))==single(Y1+(i-1)*h)
                        CNXO(i)=r;
                        break;
                    end
                end
            end
            U_enter=U_nd(:,CNXI(1:N_H));
            U_exit=U_nd(:,CNXO(1:N_H));
            U_mid=(U_enter(1,:)+U_exit(1,:))/2;
            U_max_mid=max(U_mid(1,:));
            for i=1:N_H
                U_ana_under(i,1)=sym_para(Y1,Y2,U_max_mid,Y(CNXO(i)));
            end
            figure;
            plot(U_ana_under(1:N_H),Y(CNXO(1:N_H)),U_mid(1,1:N_H),Y(CNXO(1:N_H)));
            % Pressure
            CNYC=0;
            CYC=(Y2+Y1)/2;
            for i=1:N_L
                for r=1:N
                    if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CYC)
                        CNYC(i)=r;
                        break;
                    end
                end
            end
            if CNYC==0 % Not found
                for i=1:N_L
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CYC-h/2)
                            CNYC1(i)=r;
                            break;
                        end
                    end
                end
                for i=1:N_L
                    for r=1:N
                        if single(X(r))==single(X1+(i-1)*h) && single(Y(r))==single(CYC+h/2)
                            CNYC2(i)=r;
                            break;
                        end
                    end
                end
            end
            if CNYC==0
                Rho_l=(Rho_nd(1,CNYC1(1))+Rho_nd(1,CNYC2(1)))/2;
                Rho_r=(Rho_nd(1,CNYC1(end))+Rho_nd(1,CNYC2(end)))/2;
                a=(Rho_l-Rho_r)/(X1-X2);
                b=Rho_l-a*X1;
                for i=1:N_L
                    CY(i)=X(CNYC1(i));
                    Rho_ana(i)=a*CY(i)+b;
                    Rho_d(i)=(Rho_nd(1,CNYC1(i))+Rho_nd(1,CNYC2(i)))/2;
                end
            else
                Rho_l=Rho_nd(1,CNYC(1));
                Rho_r=Rho_nd(1,CNYC(end));
                a=(Rho_l-Rho_r)/(X1-X2);
                b=Rho_l-a*X1;
                for i=1:N_L
                    CY(i)=X(CNYC(i));
                    Rho_ana(i)=a*CY(i)+b;
                    Rho_d(i)=Rho_nd(1,CNYC(i));
                end
            end
            figure;
            plot(CY,Rho_ana,CY,Rho_d);
            %%%%%%%% Error
            %%% U_y
            %%% outlet
            figure(21);
            for i=1:N_H
                U_err_y_out(i)=U_nd(2,CNXO(i));
            end
            AnaCout=zeros(1,N_H);
            plot(U_err_y_out(1:N_H),Y(CNXO(1:N_H)),AnaCout,Y(CNXO(1:N_H)));
            %%% U_x
            %%% Outlet
%             U_x_avg_in=0;       
%             for i=1:N_H-1
%                     U_x_avg_in=U_x_avg_in+(U_nd(1,CNXI(i))+U_nd(1,CNXI(i+1)))/2*h;
%             end
%             U_x_avg_in=U_x_avg_in/(Y2-Y1);
            % 2099
            U_x_avg_in=U_in(1,4)/3*2;
            U_x_avg_out=U_x_avg_in*Rho_l/Rho_r;
            U_ana_max_out=U_x_avg_out/2*3;
            for i=1:N_H
                U_ana_out_over(i,1)=sym_para(Y1,Y2,U_ana_max_out,Y(CNXO(i)));
            end
            for i=1:N_H
                U_x_out(i)=U_nd(1,CNXO(i));
                U_err_x_out(i)=U_ana_out_over(i,1)-U_x_out(i);
            end
            figure(22);
            plot(U_err_x_out(1:N_H)/U_ana_max_out,Y(CNXO(1:N_H)));
            figure(23);
            plot(U_ana_out_over(1:N_H),Y(CNXO(1:N_H))/(Y2-Y1),U_x_out(1:N_H),Y(CNXO(1:N_H))/(Y2-Y1));
            %%% Rho
            for i=1:N_L
                Rho_err(i)=(Rho_ana(1,i)-Rho_d(1,i));
            end
            figure(24);
            plot(CY,Rho_err/(abs(Rho_l-Rho_r)));
            %%% Temperature
            for i=1:N_H
                CX=(X2-X1)/2;
                if mod(N_H,2)==1
                    for r=1:N
                        if single(X(r))==single(CX) && single(Y(r))==single(Y1+(i-1)*h)
                            CNX(i)=r;
                            break;
                        end
                    end
                    Temperature=zeros(N_H,1);
                    y=zeros(N_H,1);
                else
                    for r=1:N
                        if single(X(r))==single(CX) && single(Y(r))==single(Y1+h/2+(i-1)*h)
                            CNX(i)=r;
                            break;
                        end
                    end
                    Temperature=zeros(N_H+1,1);
                    y=zeros(N_H+1,1);
                end
            end

            for i=1:N_H+1
                if i==1
                    Temperature(i)=T_bc(3);
                    y(i,1)=0;
                elseif i<N_H+1
                    Temperature(i)=T_nd(CNX(i-1));
                    y(i)=Y(CNX(i-1))/(Y2-Y1);
                else
                    Temperature(i)=T_bc(1);
                    y(i)=Y2/(Y2-Y1);
                end
            end
            
            
            T_ana_step=(T_bc(1)-T_bc(3))/(length(Temperature)-1);
            T_ana=T_bc(3):T_ana_step:T_bc(1);
            figure(25);
            plot(T_ana,y,Temperature,y);
        elseif FM==1 % arbitrary mesh
            ;
        else
            error('Wrong flag for flow type!');
        end
    elseif strcmp(top,'Moving Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1% cavity flow
        if FM==0 % IRT mesh
            U_cav=U_m(:,1);
            % U vs. Y
            figure;
            for i=1:N_H
                CX=(X2-X1)/2;
                CXL=floor(N_L/2)*h-h;
                CXR=floor(N_L/2)*h+h;
                if mod(N_H,2)==1
                    for r=1:N
                        if single(X(r))==single(CX) && single(Y(r))==single(Y1+(i-1)*h)
                            CNX(i)=r;
                            break;
                        end
                    end
                else
                    for r=1:N
                        if single(X(r))==single(CX) && single(Y(r))==single(Y1+h/2+(i-1)*h)
                            CNX(i)=r;
                            break;
                        end
                    end
                end
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
            u=zeros(N_H+1,1);
            y=zeros(N_H+1,1);
            for i=1:N_H+1
                if i==1
                    u(i)=0;
                    y(i,1)=0;
                elseif i<N_H+1
                    u(i)=U_nd(1,CNX(i-1))/U_cav(1);
                    y(i)=Y(CNX(i-1))/(Y2-Y1);
                else
                    u(i)=U_cav(1)/U_cav(1);
                    y(i)=Y2/(Y2-Y1);
                end
            end
            plot(u,y);
            % V vs. X
            figure;
            for i=1:N_L
                CY=(Y2-Y1)/2;
                if mod(N_L,2)==1
                    for r=1:N
                        if single(Y(r))==single(CY) && single(X(r))==single(X1+(i-1)*h)
                            CNY(i)=r;
                            break;
                        end
                    end
                else
                    for r=1:N
                        if single(Y(r))==single(CY) && single(X(r))==single(X1+h/2+(i-1)*h)
                            CNY(i)=r;
                            break;
                        end
                    end
                end
            end
            v=zeros(N_L+1,1);
            x=zeros(N_L+1,1);
            for i=1:N_L+1
                if i==1
                    v(i)=0;
                    x(i,1)=0;
                elseif i<N_H+1
                    v(i)=U_nd(2,CNY(i-1))/U_cav(1);
                    x(i)=X(CNY(i-1))/(X2-X1);
                else
                    v(i)=0;
                    x(i)=X2/(X2-X1);
                end
            end
            plot(x,v);
        elseif FM==1 % Arbitrary mesh
            ;
        else
            error('Wrong flag for mesh type!');
        end
    elseif strcmp(right,'Periodic')~=1 && strcmp(left,'Periodic')~=1% Other channel flow
        if N_I_N~=0 % There is immersed boundary
            Angle=zeros(1,N_I_N+1);
            angle_temp=zeros(1,N_I_N+1);
            p_cir=zeros(1,N_I_N+1);
            p_cir_temp=zeros(1,N_I_N+1);
            for r=1:N_I_N
                ND=NODE{r};
                if ND{2}>0
                    error('The node is not on the immersed boundary!');
                end
                %                 dis_horizontal=(X(r)-3)/0.3;
                dis_horizontal=(X(r)-(X2-X1)/2)/0.3;
                dis_vertical=(Y(r)-(Y2-Y1)/2)/0.3;
                angle1=acos(dis_horizontal)/pi*180;
                angle2=asin(dis_vertical)/pi*180;
                if angle2>=0 % from 0 to 180
                    if angle1<=90 && angle2<=90 % from 0 to 90
                        Angle(r)=angle2;
                    else % from 90 to 180
                        Angle(r)=angle1;
                    end
                else % from 180 to 360
                    if angle1>=90 % from 180 to 270
                        Angle(r)=180+abs(angle2);
                    else % from 270 to 360
                        Angle(r)=360-abs(angle2);
                    end
                end
                q=length(V(1,:));
                if q==7
                    p_cir(r)=Rho_nd(r)/4;
                elseif q==9
                    p_cir(r)=Rho_nd(r)/3;
                elseif q==13
                    p_cir(r)=Rho_nd(r)/2;
                else
                    error('Other lattice is not available!');
                end
            end
            front_stag_found=0;
            for s=1:r
                if single(Angle(s))==single(180)
                    front_stag_found=1;
                    break;
                end
            end
            if front_stag_found
                % Moving the front stagnation point to the first position of the
                % data and rearrange all data input
                angle_temp(1:r-s+1)=Angle(s:r)-180;
                angle_temp(r-s+2:r)=Angle(1:s-1)+180;
                p_cir_temp(1:r-s+1)=p_cir(s:r);
                p_cir_temp(r-s+2:r)=p_cir(1:s-1);
                Angle=angle_temp;
                p_cir=p_cir_temp;
            end
            % filling the last position
            Angle(r+1)=360+Angle(1);
            p_cir(r+1)=p_cir(1);
            figure
            plot(Angle,p_cir);
        end
    else
        msg=['The current type of flow cannot be plotted!'];
        disp(msg);
    end
elseif FF==1 % Taylor vortex flow
%     figure;
%     for k=1:2*a  % 2*a for showing one point; 4*a for showing 2 points
%         plot(TT,log10(U_NDY(k,:)));
%         hold on
%         plot(TT,log10(U_TV_NDY(k,:)));
%         hold on
%     end
%     hold off
end

% figure(23);
% for k=1:a
%     plot(TT,(abs(U_NDY(1+4*(k-1),:))-abs(U_NDY(3+4*(k-1),:))),TT,(abs(U_NDY(2+4*(k-1),:))-abs(U_NDY(4+4*(k-1),:))));
%     hold on
% end
% hold off
% 
% if FS==0
%     figure(24);
%     plot(TT,FCOL_c(1,:)+FCOL_c(2,:)+FCOL_c(3,:)+FCOL_c(4,:)+FCOL_c(5,:)+FCOL_c(6,:)+FCOL_c(7,:)+FCOL_c(8,:)+FCOL_c(9,:));
% end
%%%%%%%%%%%%%% Colourful vector plot
% U_min=min(U_M(1:M));
% U_rr=max(U_M(1:M))-U_min;
% VEC_kid = get(VEC,'children');
% X_kid = get(VEC_kid(1),'XData');
% Y_kid = get(VEC_kid(1),'YData');
% axis off;
% 
% figure('Position',[10 10 1000 600],'Color','w');
% cmap = jet(100); %colormap
% 
% for ii = 1:3:length(X_kid)-1
% 
%     U_c = floor((U_M(floor(ii/3)+1)-U_min)/U_rr*99)+1; %get the angle
%     ah = annotation('arrow',...
%         'Color', cmap(U_c,:),...
%         'headStyle','cback1','HeadLength',U_c/2,'HeadWidth',5);
%     set(ah,'parent',gca);
%     set(ah,'position',[X_kid(ii) Y_kid(ii) X_kid(ii+1)-X_kid(ii) Y_kid(ii+1)-Y_kid(ii)]);
% end
% axis equal tight
% axis off;