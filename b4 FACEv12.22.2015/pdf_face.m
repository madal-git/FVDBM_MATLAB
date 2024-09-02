function f_f = pdf_face(fs_dd,fs_d,fs_u,fs_uu,s,marker,dt,f_normal_latt_dot,FX,FTVD)
% function f_f = pdf_face(pdf_stencil,s,marker,FX) calculates the pdf
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
% 1---Lax-Wendroff
% 2---2nd-order upwind (SOU)
% 3---TVD
% 4---QUICK
% 5---QUICKEST
% FTVD is the flag to select a TVD scheme


% marker=upwind_marker;
if FX==0 % 1st-order upwind (SOU)
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_f=f_u;
elseif FX==1 % Lax-Wendroff
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    a=s{7}; % This alias will speed up the code
    b=s{8}; % This alias will speed up the code
%     step1=(abs(f_normal_latt_dot).*a)*dt/2;
%     step2=(b-step1)';
%     step3=step2.*(f_d-f_u);
%     f_f=f_u+step3;
    
    f_f=f_u+(b-(abs(f_normal_latt_dot).*a)*dt/2)'.*(f_d-f_u);
elseif FX==2 % 2nd-order upwind (SOU)
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{10};
    b=s{11};
    f_f=a'.*f_u+b'.*f_uu;
elseif FX==3 % TVD
    ;
elseif FX==4 % QUICK
    f_d=fs_d.*marker+fs_u.*abs(marker-1);
    f_u=fs_d.*abs(marker-1)+fs_u.*marker;
    f_uu=fs_dd.*abs(marker-1)+fs_uu.*marker;
    a=s{12};
    b=s{13};
    c=s{14};
    f_f=a'.*f_d+b'.*f_u+c'.*f_uu;
elseif FX==5 % QUICKEST
    error('Temporarily not available!');
else
    error('Wrong flag for flux calculation!');
end