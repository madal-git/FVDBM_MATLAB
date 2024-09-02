function m_eq = mrt_2d_eqm_moment (d,q,Rho,U)
%% m_eq = mrt_2d_eqm_moment (d,q,Rho,U) generetes the equilibrium moment vector based on its definition
%  d is the number of dimensions
%  q is the number of velocity components in the velocity space
%  Rho is the density
%  U is the velocity
%  m_eq is the the equilibrium moment vector
if length(Rho)~=1
    error('Check the density input');
end
if d==1
    error('Temporarily not available!');
elseif d==2
    if length(U)~=2
        error('Check the velocity input');
    end
    
    if q==5
        error('Temporarily not available!');
    elseif q==7
        error('Temporarily not available!');
    elseif q==9
        m_eq=Rho*[1, -2+3*(U'*U), 1-3*(U'*U), U(1,1), -U(1,1), U(2,1), -U(2,1), U(1,1)^2-U(2,1)^2, U(1,1)*U(2,1)]';
    else
        error('Wrong lattice model');
    end
elseif d==3
    error('Temporarily not available!');
else
    error('Wrong dimensions!');
end