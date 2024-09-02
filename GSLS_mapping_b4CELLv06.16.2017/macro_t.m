function T=macro_t(g,V,Rho)
% T=macro_t(f,V,Rho) calculates the macro temperature for Passive-Scalar model at
% each point of entire domain.
% g is the thermal pdf matrix of the entire domain
% V is the currently used lattice for thermal model
% Rho is the density vector of the entire domain

% T is the yielt temperature vector of the entire domain
R=8.314; % Thermal constant


d=length(V(:,1));
q=length(V(1,:));
p=length(g(:,1));
m=length(g(1,:));
if d~=2
    error('The current model could only perform 2-D simulations!');
end
if p~=q
    error('The lattice does not match the pdf dimension!');
end

%%%%Empty Temperature at cell centriods and nodes
T=ones(1,m);
c=ones(1,q);
for r=1:m
%     T(1,r)=c*g(:,r)/R/Rho(1,r);
    T(1,r)=c*g(:,r)/Rho(1,r);
end