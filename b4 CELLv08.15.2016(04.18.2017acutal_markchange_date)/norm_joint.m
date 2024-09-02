function C_j=norm_joint(C_b,V,C_o)
% function C_j=norm_joint(C_b,V,C_o) returns x and y coordinates of a point
% that is the joint point of two lines. The first line has direction
% defined by V and crosses the point C_b; the second line crosses C_o and
% is perpendicular to the first line.


% %% Algorithm 1, bisecing
% Tol=1e-6;
% 
% V=V/sqrt(V(1,1)^2+V(2,1)^2);
% n=[C_o(1,1)-C_b(1,1);C_o(2,1)-C_b(2,1)]; % a vector pointing to C_o from C_b, position cannot switched.
% %% make sure n and V form a angle smaller than 90 degrees
% if n'*V>=0
%     ;
% else
%     V=-V;
% end
% 
% d=dis(C_b,C_o); % distance between C_b and C_o
% C_b_t=C_b; % The other point in line defined by V
% n_t=[C_o(1,1)-C_b_t(1,1);C_o(2,1)-C_b_t(2,1)]; % Form the third line
% while n_t'*V>=0
%     C_b_t=C_b_t+d*V; % Along line one to find the next C_b_t that satisfy the
%     % angle between first line and third line is greater than 90 degrees
%     % then we know C_j is between C_b and C_b_t
%     n_t=[C_o(1,1)-C_b_t(1,1);C_o(2,1)-C_b_t(2,1)]; % Form the new third line
% end
% %% Use bisection method to find C_j
% C_l=C_b;
% C_h=C_b_t;
% C_m=(C_l+C_h)/2;
% n_j=[C_o(1,1)-C_m(1,1);C_o(2,1)-C_m(2,1)];
% while abs(n_j'*V)>=Tol
%     if n_j'*V>0
%         C_l=C_m;
%     else
%         C_h=C_m;
%     end
%     C_m=(C_l+C_h)/2;
%     n_j=[C_o(1,1)-C_m(1,1);C_o(2,1)-C_m(2,1)];
% end
% %% Return
% C_j=C_m;

%% Algorithm 2, using function intercept
Tol2=1e-8;
V=V/sqrt(V(1,1)^2+V(2,1)^2);
L=dis(C_b,C_o);
C1=C_b+10*L*V;
C2=C_b-10*L*V;
if abs(V(2,1))<Tol2
    m1=0;
    m2=1;
else
    m1=1;
    m2=-V(1,1)/V(2,1)*m1;
end
V_n=[m1;m2];
[Found, C_i]=intercept(C_o,V_n,C1,C2);
if Found
    C_j=C_i;
else
    error('Logic error!');
end
