function F = force(lattice, weight, Density, Temperature, T_ref, gravity, FFC, FTH)
%% function F = force(lattice, weight, Density, Temperature, T_ref, gravity)
%  calculates the vector of force in each lattice direction for all centroids from known conditions

% F is the resulted force vector in each direction of the selected lattice
% lattice is the selected lattice
% weight is the corresponding weights
% Density is the density matrix of the entire domain
% Temperature is the temperature of the entire domain
% T_ref is the reference temperature, a single value
% gravity is the gravititional force, which is a sigle vector
% FFC is the flag for the force term flag FFC=0---No force term; FFC=1---Force term
% FTH is the flag for thermal problem. FTH=0---isothermal; FTH=1---Thermal

%% Checking
q1=length(lattice);
q2=length(weight);
if q1~=q2
    error('The selected lattice and its weight donnot match!');
else
    q=q1;
end

M1=length(Density);
M2=length(Temperature);
if M1~=M2
    error('The size of the computational domain for density and temperature donnot match!');
else
    M=M1;
end

%% Determing speed of sound
if q==9
    c_s=1/sqrt(3);
else
    error('Other lattices are not available!');
end

%% Calculating force
if FFC==0
    error('Wrong flag for force term!');
else
    F_source = body_force (gravity, Density, Temperature, T_ref, FTH, 1);
end

w=zeros(q,M);
for i=1:M
    w(:,i)=weight';
end

F=((lattice')*F_source).*w/c_s^2;



