function Phi=TVD(r,R,Ftvd)
% function Phi=TVD(r,R,Ftvd) generates the Phi function
% based on variety of TVD scheme to calculate the 2nd-order upwind flux
% on non-uniform three-point stencil
% r is the ratio of gradient from upwind point to further upwind point 
% versus the gradient from downwind point to upwind point
% R is the spatial ratio. R=(downwind cell size + upwind cell size)/upwind cell size
% Ftvd is the flag for different TVD scheme
L1=length(r);
L2=length(R);
if L1~=L2
    error('r and R should have the same dimensions!');
end
if Ftvd==0 % minmod
    Phi=max(zeros(L1,1),min(r,2*ones(L1,1)));
elseif Ftvd==1 % SUPERBEE
    Phi_temp=max(min(R.*r,ones(L1,1)),min(r,R));
    Phi=max(zeros(L1,1),Phi_temp);
elseif Ftvd==2 % Van Leer
    Phi=zeros(L1,1);
    for i=1:L1
        if r(i,1)==Inf
            Phi(i,1)=R(i,1);
        else
            Phi(i,1)=(R(i,1)*r(i,1)/2+R(i,1)*abs(r(i,1))/2)/(R(i,1)+r(i,1)-1);
        end
    end
elseif Ftvd==3 % Osher
    Phi=max(zeros(L1,1),min(r,R));
elseif Ftvd==4 % Lax-Wendroff or central difference
    Phi=ones(L1,1);
else
    error('The flag for TVD scheme is unavailable or invalid!');
end


% Method 1
% for i=1:L1
%     if Phi(i,1)<0
%         Phi(i,1)=0;
%     end
% end
% Method 2
Phi=Phi.*(Phi>=0);