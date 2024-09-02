function [RHO,U,T]=macro_h(f,V,Rho_r,fd)
% [RHO,U]=macro_hydro(f,V) calculates the macro density and velocity at
% each point of entire domain.
% f is the pdf matrix of the entire domain
% V is the currently used lattice
% Rho_f is the reference density
% fd is the flag for which density is used. fd=0----local density; fd=1----reference density
% RHO is yielt density vector of the entire domain
% U is the yielt velocity matrix of the entire domain
% T is the yielt temperature vector of the entire domain

d=length(V(:,1));
q=length(V(1,:));
p=length(f(:,1));
m=length(f(1,:));
if d~=2
    error('The current model could only perform 2-D simulations!');
end
if p~=q
    error('The lattice does not match the pdf dimension!');
end
%%%%Empty density at cell centriods and nodes
RHO=zeros(1,m);
%%%%Empty velocity at cell centriods and nodes
U=zeros(d,m);
%%%%Empty Temperature at cell centriods and nodes
T=ones(1,m);
for r=1:m
    RHO(1,r)=sum(f(:,r));
    U(:,r)=V*f(:,r);
    if fd==0
        U(:,r)=U(:,r)/RHO(1,r);
    elseif fd==1
        U(:,r)=U(:,r)/Rho_r;
    else
        error('Flag for density is invalid or unavailable!');
    end
    %% Temperature
    if q==37
        c=V-U(:,r)*ones(1,q);
        T(1,r)=f(:,r)'*diag(c'*c)/2/RHO(1,r);
    end
end