function f_c = pdf_face_new(fs_dd,fs_d,fs_u,fs_uu,fs_dd_ctd,fs_d_ctd,fs_u_ctd,fs_uu_ctd,s,marker,dt,f_normal_latt_dot,f_c_old,FX,FUPD,FTVD)
% function f_c = pdf_face(pdf_stencil,s,marker,FX) calculates the pdf
% values at the face center depending on the pdf values along the stencil
% points of the face and the given flag for certian scheme of flux
% calculation
% f_f is the flux on the face
% pdf_stencil is the pdf values at each stencil point [dd,d,u,uu] (based on the face normal)
% s is the assorted stencil info characterized by the applied lattice
% structure
% marker is n_by_1 vector to tell the direction of each lattice velocity
% with respect to the face normal. 1---dot(face normal, lattice velocity)>=0;
%                                  0---dot(face normal, lattice velocity)<0
% dt is the time step size
% f_normal_latt_dot is the dot product of face normal and lattice velocities
% FX is the flag to select a scheme of flux calculation
% 0---1st-order upwind (FOU)
% 1---Godunov-PL, Lax-Wendroff
% 2---Godunov-PL, Beam-Warming
% 3---Godunov-PL, Fromm
% 4---Non-Godunov general form
% 5---Non-Godunov,2nd-order upwind (SOU)
% 6---Non-Godunov, TVD
% 7---Non-Godunov,QUICK
% 8---Non-Godunov,QUICKEST
% 9---Godunov-PP
% 10--Godunov, TVD, PC + PL-LW
% 11--Godunov, TVD, PC + PP
% 12--Godunov, TVD, PL-LW + PP
% FUPD is the flag for wether use upwind scheme to calculate the flux
% FTVD is the flag to select a TVD scheme
% marker=upwind_marker;

e=10;
if FX==0 % 1st-order upwind (FOU) or PC Godunov
    if FUPD==0 % No upwind
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
        f_c=f_u;
    elseif FUPD==1 % upwind
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
        f_c=f_u;
    else
        error('Wrong flag for FUPD');
    end
    % no difference between upwind and no upwind
elseif FX==1 % Lax-Wendroff, PL Godunov
    if FUPD==0 % No upwind
%         f_d=fs_d.*marker+fs_u.*abs(marker-1);
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
        one_over_L_u=s{7}./s{8};
        %     step1=(abs(f_normal_latt_dot).*a)*dt/2;
        %     step2=(b-step1)';
        %     step3=step2.*(f_d-f_u);
        %     f_f=f_u+step3;
        
%         f_c=f_u+(1-(abs(f_normal_latt_dot).*one_over_L_u)*dt/2)'.*(f_c_old-f_u); %% Too diffusive
        f_c=(f_u+(1-(abs(f_normal_latt_dot).*one_over_L_u)*dt/2)'.*(f_c_old/2-f_u))./(0.5+(abs(f_normal_latt_dot).*one_over_L_u)*dt/4)';
    elseif FUPD==1 % upwind
        f_d=fs_d_ctd.*marker+fs_u_ctd.*abs(marker-1);
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%         f_d=fs_d.*marker+fs_u.*abs(marker-1);
%         f_u=fs_d_ctd.*abs(marker-1)+fs_u_ctd.*marker;
        a=s{7}; % This alias will speed up the code
        b=s{8}; % This alias will speed up the code
        %     step1=(abs(f_normal_latt_dot).*a)*dt/2;
        %     step2=(b-step1)';
        %     step3=step2.*(f_d-f_u);
        %     f_f=f_u+step3;
        f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_d-f_u);
    else
        error('Wrong flag for FUPD');
    end
elseif FX==2 % Beam-Warming
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{9};
    b=s{10};
    f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_u-f_uu);
elseif FX==3 % Fromm
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{11};
    b=s{12};
    f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_d-f_uu);
elseif FX==4 % Face-value-based scheme
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    FFC=1;
    if FFC==1 % Central difference
        c1=s{4}.*s{7};
        c2=s{8};
        f_c_plus=c1'.*f_u+c2'.*f_d;
    elseif FFC==2 % SOU
        c1=s{16};
        c2=s{17};
        f_c_plus=c1.*f_u+c2.*f_uu;
    elseif FFC==3 % QUICK
        c1=s{18};
        c2=s{19};
        c3=s{20};
        f_c_plus=c1.*f_d+c2.*f_u+c3.*f_uu;
    else
        error('Wrong flag!');
    end
    a=s{13};
    b=s{14};
    f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_c_plus-f_u);
elseif FX==5 % 2nd-order upwind (SOU)
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{16};
    b=s{17};
    f_c=a.*f_u+b.*f_uu;
elseif FX==6 % TVD, non-Godunov
    e=10;
    % Method 1,vectorization
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    % Method 2, for loop
%     L=length(marker);
%     f_d=zeros(L,1);
%     f_u=f_d;
%     f_uu=f_d;
%     for i=1:L
%         if marker(i,1)==1
%             f_d(i,1)=fs_d(i,1);
%             f_u(i,1)=fs_u(i,1);
%             f_uu(i,1)=fs_uu(i,1);
%         else
%             f_d(i,1)=fs_u(i,1);
%             f_u(i,1)=fs_d(i,1);
%             f_uu(i,1)=fs_dd(i,1);
%         end
%     end
    
    a=s{8};
    b=s{15};
    L=length(f_d);
    %% r parameter
    r=zeros(L,1);
    for i=1:L
        if double(e+f_d(i,1))==double(e+f_u(i,1))
            if double(e+f_u(i,1))==double(e+f_uu(i,1))
                r(i,1)=1;
            else
                r(i,1)=Inf;
            end
        else
            r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
        end
    end
    f_c=f_u+((f_d-f_u).*TVD(r,(1./a)',FTVD)).*a';
elseif FX==7 % QUICK
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{18};
    b=s{19};
    c=s{20};
    f_c=a.*f_d+b.*f_u+c.*f_uu;
elseif FX==8 % QUICKEST
    error('Temporarily not available!');
elseif FX==9 % Third-order Godunov, PP
    %% Algorithm 1, use f_d, f_u and f_uu, even more diffusion error than the linear reconstruction
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=s{24};
%     A12=s{21};
%     A13=ones(L,1);
%     
%     A21=s{25};
%     A22=s{22};
%     A23=A13;
%     
%     A31=s{26};
%     A32=s{23};
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_d(i,1);f_u(i,1);f_uu(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3-b.*abs(f_normal_latt_dot)'*dt/2+c;
    %% Algorithm 2, use f_c+, f_u and f_uu, more diffusion error than linear reconstruction
    % Calculate the coefficients for parabolic reconstruction
    if FUPD==0 % No upwind
        f_d=fs_d.*marker+fs_u.*abs(marker-1);
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
        f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
        f_c=s{28}.*f_u+s{29}.*f_d;
        
        %     % Method 1, do matrix inverse with the magic backslash
        %     L=length(f_d);
        %     a=zeros(L,1);
        %     b=a;
        %     c=a;
        %     A11=zeros(L,1);
        %     A12=zeros(L,1);
        %     A13=ones(L,1);
        %
        %     A21=s{25};
        %     A22=s{22};
        %     A23=A13;
        %
        %     A31=s{26};
        %     A32=s{23};
        %     A33=A13;
        %     for i=1:L
        %         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c(i,1);f_u(i,1);f_uu(i,1)];
        %         a(i,1)=C(1,1);
        %         b(i,1)=C(2,1);
        %         c(i,1)=C(3,1);
        %     end
        % Method 2, get the analytic solution of the matrix inverse
        a=s{21}.*f_c+s{22}.*f_u+s{23}.*f_uu;
        b=s{24}.*f_c+s{25}.*f_u+s{26}.*f_uu;
        c=f_c;
        % Calculate face value
        f_c=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    elseif FUPD==1 % Upwind
        f_d=fs_d_ctd.*marker+fs_u_ctd.*abs(marker-1);
        f_u=fs_d.*abs(marker-1)+fs_u.*marker;
        f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
        f_c=s{28}.*f_u+s{29}.*f_d;
        a=s{21}.*f_c+s{22}.*f_u+s{23}.*f_uu;
        b=s{24}.*f_c+s{25}.*f_u+s{26}.*f_uu;
        c=f_c;
        % Calculate face value
        f_c=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    else
        error('Wrong flag for FUPD');
    end
    
    %% Algorithm 3, use f_c+, f_c- and f_u, more diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus=(f_u+f_uu)/2;
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 4, use f_c+, f_c- and f_u, where f_c- takes the value of f_uu but keep the location
    % it is supposed to be, less diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus=f_uu;
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 5, use f_c+, f_c- and f_u, where f_c- takes the value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, less diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 6, use f_c+, f_c- and f_u, where f_c- takes the second value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, more diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus1-(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 7, use f_c+, f_c- and f_u, where f_c- takes the third type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, less diffusive than Algorithm 2, so far the best
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+2*(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 8, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be. It is the genral form of Algorithm 5 and 7
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     
% 
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2*FPPI-f_c_minus1*(FPPI-1);
% %     f_c_minus=f_c_minus2*(2-FPPI)+f_c_minus1*(FPPI-1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 9, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 with judgement of r parameter but keep the location
%     % it is supposed to be.
%     % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     % Calculate the r parameter
%     b=s{15};
%     L=length(f_d);
%     r=zeros(L,1);
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
%         end
%     end
%     % Calculate the guessed values at positive and negative cell interface
%     % positive interface is linearly interpolated
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     % negative interface is linealy interpolated and then tuned based on r
%     % parameter
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=zeros(L,1);
%     for i=1:L
%         if r(i,1)<0  % A=1 for f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1); A=2 for f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1)
%             f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
% %             f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%         else
%             if r(i,1)<=1  % B
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%             else   % C
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%             end
%             
%         end
%     end
%     % PP reconstruction
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 10, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 with judgement of r parameter but keep the location
    % it is supposed to be. It replace the FPPI with the function of r parameter
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     % Calculate the r parameter
%     b=s{15};
%     L=length(f_d);
%     r=zeros(L,1);
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
%         end
%     end
%     R=abs(r);
% %     R=1./abs(r);
%     for i=1:L
%         if R(i,1)>2
%             R(i,1)=2;
%         end
%     end
%     % Calculate the guessed values at positive and negative cell interface
%     % positive interface is linearly interpolated
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     % negative interface is linealy interpolated and then tuned based on r
%     % parameter
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=zeros(L,1);
%     for i=1:L
%         if r(i,1)<0  % A=1 for f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1); A=2 for f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1)
%             f_c_minus(i,1)=f_c_minus2(i,1)*R(i,1)-f_c_minus1(i,1)*(R(i,1)-1);
% %             f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%         else
%             if r(i,1)<=1  % B
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-R(i,1))+f_c_minus1(i,1)*(R(i,1)-1);
%             else   % C
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-R(i,1))+f_c_minus1(i,1)*(R(i,1)-1);
%             end
%             
%         end
%     end
%     % PP reconstruction
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 11, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 with judgement of r parameter but keep the location
    % it is supposed to be. It replace the FPPI with the function of r parameter
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     % Calculate the r parameter
%     b=s{15};
%     L=length(f_d);
%     r=zeros(L,1);
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
%         end
%     end
%     % Calculate the guessed values at positive and negative cell interface
%     % positive interface is linearly interpolated
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     % negative interface is linealy interpolated and then tuned based on r
%     % parameter
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=zeros(L,1);
%     for i=1:L
%         if r(i,1)==1  % A=1 for f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1); A=2 for f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1)
%             f_c_minus(i,1)=f_c_minus1(i,1);
% %             f_c_minus(i,1)=f_c_minus2(i,1);
%         else % B
%             f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%         end
%     end
%     % PP reconstruction
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 12, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 with judgement of r parameter but keep the location
    % it is supposed to be. It replace the FPPI with the function of r parameter
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     % Calculate the r parameter
%     b=s{15};
%     L=length(f_d);
%     r=zeros(L,1);
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
%         end
%     end
% %     R=abs(r);
% % %     R=1./abs(r);
% %     for i=1:L
% %         if R(i,1)>2
% %             R(i,1)=2;
% %         end
% %     end
%     % Calculate the guessed values at positive and negative cell interface
%     % positive interface is linearly interpolated
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     % negative interface is linealy interpolated and then tuned based on r
%     % parameter
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=zeros(L,1);
%     for i=1:L
%         if r(i,1)<0  % A=1 for f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1); A=2 for f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1)
% %             f_c_minus(i,1)=f_c_minus2(i,1)*R(i,1)-f_c_minus1(i,1)*(R(i,1)-1);
%             f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%         else
%             if r(i,1)==1  % B
%                 f_c_minus(i,1)=f_c_minus1(i,1);
%             else   % C
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-R(i,1))+f_c_minus1(i,1)*(R(i,1)-1);
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%             end
%             
%         end
%     end
%     % PP reconstruction
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 13, use f_c+, f_c- and f_u, where f_c- takes the general form of blended f_uu and (f_uu+f_u)/2 with judgement of r parameter but keep the location
    % it is supposed to be. It replace the FPPI with the function of r parameter
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     % Calculate the r parameter
%     b=s{15};
%     L=length(f_d);
%     r=zeros(L,1);
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*b(i,1);
%         end
%     end
%     R=abs(r);
% %     R=1./abs(r);
%     for i=1:L
%         if R(i,1)>2
%             R(i,1)=2;
%         end
%     end
%     % Calculate the guessed values at positive and negative cell interface
%     % positive interface is linearly interpolated
%     f_c_plus=s{28}.*f_u+s{29}.*f_d;
%     % negative interface is linealy interpolated and then tuned based on r
%     % parameter
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=zeros(L,1);
%     for i=1:L
%         if r(i,1)<0  % A=1 for f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1); A=2 for f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1)
%             f_c_minus(i,1)=f_c_minus2(i,1)*R(i,1)-f_c_minus1(i,1)*(R(i,1)-1);
% %             f_c_minus(i,1)=f_c_minus2(i,1)*FPPI-f_c_minus1(i,1)*(FPPI-1);
%         else
%             if r(i,1)==1  % B
%                 f_c_minus(i,1)=f_c_minus1(i,1);
%             else   % C
%                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-R(i,1))+f_c_minus1(i,1)*(R(i,1)-1);
% %                 f_c_minus(i,1)=f_c_minus2(i,1)*(2-FPPI)+f_c_minus1(i,1)*(FPPI-1);
%             end
%             
%         end
%     end
%     % PP reconstruction
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 14, use f_c+, f_c- and f_u, where f_c- takes the first type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, and f_c+ takes the first type value of blended f_d and (f_d+f_u)/2 but keep the location
    % it is supposed to be,
    % f_c_plus2=f_d
    % INCORRECT, UNPHYSICAL
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus1=s{28}.*f_u+s{29}.*f_d;
%     f_c_plus2=f_d;
%     f_c_plus=f_c_plus2+(f_c_plus2-f_c_plus1);
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 15, use f_c+, f_c- and f_u, where f_c- takes the first type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, and f_c+ takes the first type value of blended f_d and (f_d+f_u)/2 but keep the location
    % it is supposed to be,
    % f_c_plus2=s{28}.*f_u+s{29}.*f_d, this is actually defining f_c+=f_u
    % more diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus2=s{28}.*f_u+s{29}.*f_d;
%     f_c_plus1=f_d;
%     f_c_plus=f_c_plus2+(f_c_plus2-f_c_plus1);
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 16, use f_c+, f_c- and f_u, where f_c- takes the first type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, and f_c+ takes the first type value of blended f_u and (f_d+f_u)/2 but keep the location
    % it is supposed to be,
    % f_c_plus2=f_u
    % too much more diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus1=s{28}.*f_u+s{29}.*f_d;
%     f_c_plus2=f_u;
%     f_c_plus=f_c_plus2+(f_c_plus2-f_c_plus1);
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 17, use f_c+, f_c- and f_u, where f_c- takes the first type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, and f_c+ takes f_d but keep the location
    % it is supposed to be,
    % f_c_plus2=f_u
    % too much negative diffusivity, explodes
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus=f_d;
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 18, use f_c+, f_c- and f_u, where f_c- takes the third type value of blended f_uu and (f_uu+f_u)/2 but keep the location
    % it is supposed to be, and f_c+ takes the third type value of blended f_u and (f_d+f_u)/2 but keep the location
    % it is supposed to beless diffusive than Algorithm 2
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     f_c_plus1=s{28}.*f_u+s{29}.*f_d;
%     f_c_plus2=f_d;
%     f_c_plus=f_c_plus2+2*(f_c_plus2-f_c_plus1);
%     f_c_minus1=(f_u+f_uu)/2;
%     f_c_minus2=f_uu;
%     f_c_minus=f_c_minus2+2*(f_c_minus2-f_c_minus1);
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
    %% Algorithm 19, use f_c+, f_c- and f_u, where f_c- takes f_uu, f_c+ takes f_d, which squeeze the profile to make it steeper, less diffusive, but oscillate
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
% %     f_c_plus=s{28}.*f_u+s{29}.*f_d;
% f_c_plus=f_d;
%     f_c_minus=(f_u+f_uu)/2;
% % f_c_minus=f_uu;
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=zeros(L,1);
%     A12=zeros(L,1);
%     A13=ones(L,1);
%     
%     A21=(s{5}.*s{5})';
%     A22=-s{5}';
%     A23=A13;
%     
%     A32=-s{5}'-(s{6}-s{5})'/2;
%     A31=A32.*A32;
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_c_plus(i,1);f_u(i,1);f_c_minus(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate face value
%     f_f=(a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c)*1;
elseif FX==10 % PL with SOU for f_c
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    c1=s{16};
    c2=s{17};
    f_c=c1.*f_u+c2.*f_uu;
    
    a=s{13};
    b=s{14};
    f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_c-f_u);
elseif FX==11 % PL with QUICK for f_c
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    c1=s{18};
    c2=s{19};
    c3=s{20};
    f_c=c1.*f_d+c2.*f_u+c3.*f_uu;
    
    a=s{13};
    b=s{14};
    f_c=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_c-f_u);
elseif FX==12 % PP with SOU for f_c
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    c1=s{16};
    c2=s{17};
    f_c=c1.*f_u+c2.*f_uu;
    
    a=s{21}.*f_c+s{22}.*f_u+s{23}.*f_uu;
    b=s{24}.*f_c+s{25}.*f_u+s{26}.*f_uu;
    c=f_c;
    % Calculate face value
    f_c=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
elseif FX==13 % PP with QUICK for f_c
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    c1=s{18};
    c2=s{19};
    c3=s{20};
    f_c=c1.*f_d+c2.*f_u+c3.*f_uu;
    
    a=s{21}.*f_c+s{22}.*f_u+s{23}.*f_uu;
    b=s{24}.*f_c+s{25}.*f_u+s{26}.*f_uu;
    c=f_c;
    % Calculate face value
    f_c=a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3 - b.*abs(f_normal_latt_dot)'*dt/2 + c;
elseif FX==14 % Non linear methods, PC + limited PL
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{7}; % This alias will speed up the code
    b=s{8}; % This alias will speed up the code
    f_c_plus=s{28}.*f_u+s{29}.*f_d;
    % r parameter
    e=10;
    L=length(f_d);
    r=zeros(L,1);
    coe=s{27};
    for i=1:L
        if double(e+f_c_plus(i,1))==double(e+f_u(i,1))
            if double(e+f_u(i,1))==double(e+f_uu(i,1))
                r(i,1)=1;
            else
                r(i,1)=Inf;
            end
        else
            r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_c_plus(i,1)-f_u(i,1))*coe(i,1);
        end
    end
    % Calculate f_face
%     f_f=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_d-f_u).*TVD(r,(1./c1)',FTVD);
    f_c=f_u+((b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_d-f_u)).*TVD(r,2*ones(L,1),FTVD);
elseif FX==15 % Non linear methods, PC + limited PP
    %% Algorithm 1, use f_c, f_u and f_uu, PC + limited PP
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=s{24};
%     A12=s{21};
%     A13=ones(L,1);
%     
%     A21=s{25};
%     A22=s{22};
%     A23=A13;
%     
%     A31=s{26};
%     A32=s{23};
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_d(i,1);f_u(i,1);f_uu(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate limiter
%     e=10;
%     r=zeros(L,1);
%     coe=s{15};
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*coe(i,1);
%         end
%     end
%     % Calculate f_face
%     f_f=f_u+(a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3-b.*abs(f_normal_latt_dot)'*dt/2+c-f_u).*(TVD(r,2*ones(L,1),FTVD));
    %% Algorithm 2, use f_c+, f_u and f_uu, PC + limited PP
    % Calculate the coefficients for parabolic reconstruction
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    f_c_plus=s{28}.*f_u+s{29}.*f_d;
    a=s{21}.*f_c_plus+s{22}.*f_u+s{23}.*f_uu;
    b=s{24}.*f_c_plus+s{25}.*f_u+s{26}.*f_uu;
    c=f_c_plus;
    % Calculate limiter
    e=10;
    L=length(f_d);
    r=zeros(L,1);
    coe=s{27};
    for i=1:L
        if double(e+f_c_plus(i,1))==double(e+f_u(i,1))
            if double(e+f_u(i,1))==double(e+f_uu(i,1))
                r(i,1)=1;
            else
                r(i,1)=Inf;
            end
        else
            r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_c_plus(i,1)-f_u(i,1))*coe(i,1);
        end
    end
    % Calculate f_face
    f_c=f_u+(a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3-b.*abs(f_normal_latt_dot)'*dt/2+c-f_u).*(TVD(r,2*ones(L,1),FTVD));
elseif FX==16 % Non linear methods, PL(LW) + limited PP
    %% Algorithm 1, use f_c, f_u and f_uu, PL(LW) + limited PP
    % Calculate the coefficients for parabolic reconstruction
%     f_d=fs_d.*marker+fs_u.*abs(marker-1);
%     f_u=fs_d.*abs(marker-1)+fs_u.*marker;
%     f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
%     c1=s{7}; % This alias will speed up the code
%     c2=s{8}; % This alias will speed up the code
%     %     step1=(abs(f_normal_latt_dot).*a)*dt/2;
%     %     step2=(b-step1)';
%     %     step3=step2.*(f_d-f_u);
%     %     f_f=f_u+step3;
%     f_f_PL=f_u+(c2-(abs(f_normal_latt_dot).*c1)*dt/2)'.*(f_d-f_u);
%     
%     L=length(f_d);
%     a=zeros(L,1);
%     b=a;
%     c=a;
%     
%     A11=s{24};
%     A12=s{21};
%     A13=ones(L,1);
%     
%     A21=s{25};
%     A22=s{22};
%     A23=A13;
%     
%     A31=s{26};
%     A32=s{23};
%     A33=A13;
%     for i=1:L
%         C=[A11(i),A12(i),A13(i);A21(i),A22(i),A23(i);A31(i),A32(i),A33(i)]\[f_d(i,1);f_u(i,1);f_uu(i,1)];
%         a(i,1)=C(1,1);
%         b(i,1)=C(2,1);
%         c(i,1)=C(3,1);
%     end
%     % Calculate limiter
%     e=10;
%     r=zeros(L,1);
%     coe=s{15};
%     for i=1:L
%         if double(e+f_d(i,1))==double(e+f_u(i,1))
%             if double(e+f_u(i,1))==double(e+f_uu(i,1))
%                 r(i,1)=1;
%             else
%                 r(i,1)=Inf;
%             end
%         else
%             r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_d(i,1)-f_u(i,1))*coe(i,1);
%         end
%     end
%     % Calculate f_face
%     f_f=f_f_PL+(a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3-b.*abs(f_normal_latt_dot)'*dt/2+c-f_f_PL).*(TVD(r,2*ones(L,1),FTVD));
    %% Algorithm 2, use f_c+, f_u and f_uu, PL(LW) + limited PP
    % Calculate the coefficients for parabolic reconstruction
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    f_c_plus=s{28}.*f_u+s{29}.*f_d;
    c1=s{7}; % This alias will speed up the code
    c2=s{8}; % This alias will speed up the code
    f_f_PL=f_u+(c2-(abs(f_normal_latt_dot).*c1)*dt/2)'.*(f_d-f_u);
    
    a=s{21}.*f_c_plus+s{22}.*f_u+s{23}.*f_uu;
    b=s{24}.*f_c_plus+s{25}.*f_u+s{26}.*f_uu;
    c=f_c_plus;
    % Calculate limiter
    e=10;
    L=length(f_d);
    r=zeros(L,1);
    coe=s{27};
    for i=1:L
        if double(e+f_c_plus(i,1))==double(e+f_u(i,1))
            if double(e+f_u(i,1))==double(e+f_uu(i,1))
                r(i,1)=1;
            else
                r(i,1)=Inf;
            end
        else
            r(i,1)=(f_u(i,1)-f_uu(i,1))/(f_c_plus(i,1)-f_u(i,1))*coe(i,1);
        end
    end
    % Calculate f_face
    f_c=f_f_PL+(a.*(abs(f_normal_latt_dot).*abs(f_normal_latt_dot))'*dt^2/3-b.*abs(f_normal_latt_dot)'*dt/2+c-f_f_PL).*(TVD(r,2*ones(L,1),FTVD));
else
    error('Wrong flag for flux calculation!');
end