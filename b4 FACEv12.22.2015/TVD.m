function Phi=TVD(r,R,Ftvd)
% function Phi=TVD(r,R,Ftvd) generates the Phi function
% based on variety of TVD scheme to calculate the 2nd-order upwind flux
% on non-uniform three-point stencil
% r is the ratio of gradient from upwind point to further upwind point 
% versus the gradient from downwind point to upwind point
% R is the spatial ratio. R=(downwind cell size + upwind cell size)/upwind cell size
% Ftvd is the flag for different TVD scheme

if Ftvd==0 % minmod
    Phi=max(0,min(r,2));
elseif Ftvd==1 % SUPERBEE
    Phi_temp=max(min(R*r,1),min(r,R));
    Phi=max(0,Phi_temp);
elseif Ftvd==2 % Van Leer
    if r==Inf
        Phi=R;
    else
        Phi=(R*r/2+R*abs(r)/2)/(R+r-1);
    end
elseif Ftvd==3 % Osher
    Phi=max(0,min(r,R));
elseif Ftvd==4 % Lax-Wendroff or central difference
    Phi=1;
else
    error('The flag for TVD scheme is unavailable or invalid!');
end

if Phi<0
    Phi=0;
end