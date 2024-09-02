function g_eq=eqm_t(V,U,rho,T,q,w,Rho_r,fd)
% g_eq=eqm_t(V,U,rho,T,q,w,Rho_r,fd) calculates the thermodynamic equilibrium pdf 
% according to given lattice
% V is the velocity decomposition of the lattice
% U is the column vector of velocity at current point (could be centroid or node)
% rho is the density at current point (could be centroid or node)
% T is the temperature at current point (could be centroid or node)
% q is the total number of velocity components in the lattice
% w is the weighting factor for calculating equilibrium pdf accoring to 
% the lattice structure
% Rho_f is the reference density
% fd is the flag for which density is used. fd=0----local density; fd=1----reference density

R=8.314; % Thermal constant

if q==7 % D2Q7
    error('Temporarily not available!');
elseif q==9 % D2Q9
    g_eq=zeros(q,1);
    
    g_eq(1)=w(1)*rho*R*T*((U')*U);
    
    g_eq(2)=w(2)*rho*R*T*(3/2+3/2*V(:,2)'*U+9/2*(V(:,2)'*U)^2-3/2*(U'*U));
    
    g_eq(3)=w(3)*rho*R*T*(3/2+3/2*V(:,3)'*U+9/2*(V(:,3)'*U)^2-3/2*(U'*U));
    
    g_eq(4)=w(4)*rho*R*T*(3/2+3/2*V(:,4)'*U+9/2*(V(:,4)'*U)^2-3/2*(U'*U));
    
    g_eq(5)=w(5)*rho*R*T*(3/2+3/2*V(:,5)'*U+9/2*(V(:,5)'*U)^2-3/2*(U'*U));
    
    g_eq(6)=w(6)*rho*R*T*(3+6*V(:,6)'*U+9/2*(V(:,6)'*U)^2-3/2*(U'*U));
    
    g_eq(7)=w(7)*rho*R*T*(3+6*V(:,7)'*U+9/2*(V(:,7)'*U)^2-3/2*(U'*U));
    
    g_eq(8)=w(8)*rho*R*T*(3+6*V(:,8)'*U+9/2*(V(:,8)'*U)^2-3/2*(U'*U));
    
    g_eq(9)=w(9)*rho*R*T*(3+6*V(:,9)'*U+9/2*(V(:,9)'*U)^2-3/2*(U'*U));
elseif q==13 % D2Q13
    error('Temporarily not available!');
else
    error('Wrong lattice type for thermal model!');
end