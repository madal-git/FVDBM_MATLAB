function matrix = mrt_2d_relax_matrix (d,q,tau)
%% function matrix = mrt_2d_relax_matrix (d,q,tau) generetes the relaxation matrix and its inverse matrix for the MRT model
%  d is the number of dimensions
%  q is the number of velocity components in the velocity space
%  tau is the relaxation time
%  matrix is the relaxation matrix

if d==1
    error('Temporarily not available!');
elseif d==2
    if q==5
        error('Temporarily not available!');
    elseif q==7
        error('Temporarily not available!');
    elseif q==9
        s_v=1/tau;
%         s_e=1/tau; % Tuning
%         s_epsilon=1/tau; % Tuning
%         s_q=1/tau; % Tuning
        %% The following values are from Dominique d'Humieres. Multiple-relaxation-time lattice Boltzmann models in three dimensions. Phi. Tran. Royal Soc. A, 360(1792), 2002
        s_e=1.19; % Tuning
        s_epsilon=1.4; % Tuning
        s_q=1.2; % Tuning
        %% The following values are from Nemer Mahammad Master Thesis.
%         s_e=0.5; % Tuning
%         s_epsilon=0.8; % Tuning
%         s_q=8*(2-s_v)/(8-s_v); % Tuning
        %% The following values are from Nemer Mahammad Master Thesis.
%         s_e=1.63; % Tuning
%         s_epsilon=1.14; % Tuning
%         s_q=8*(2-s_v)/(8-s_v); % Tuning

        matrix=diag([0, s_e, s_epsilon, 0, s_q, 0, s_q, s_v, s_v]);
    else
        error('Wrong lattice model');
    end
elseif d==3
    error('Temporarily not available!');
else
    error('Wrong dimensions!');
end