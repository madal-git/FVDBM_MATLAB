function f_eq=eqm_h(V,U,rho,T,q,w,Rho_r,fd)
% f_eq=eqm_h(V,U,rho,q,w) calculates the hydrodynamic equilibrium pdf 
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

if q==7 % D2Q7
    if fd==0  
        % for loop
%         for l=1:q
%             f_eq(l,1)=w(l)*rho*(1+4*(V(:,l)')*U+8*((V(:,l)')*U)^2-2*(U')*U);
%         end
        % vectorized
        f_eq=diag(rho*(1+4*V'*U+8*diag((V'*U)*(V'*U)')-2*(U')*U)*w);
    elseif fd==1
        % for loop
%         for l=1:q
%             f_eq(l,1)=w(l)*(rho+Rho_r*(4*(V(:,l)')*U+8*((V(:,l)')*U)^2-2*(U')*U));
%         end
        % vectorized
        f_eq=diag(rho+Rho_r*(4*V'*U+8*diag((V'*U)*(V'*U)')-2*(U')*U)*w);
    else
        error('Flag for density is invalid or unavailable!');
    end
elseif q==9 % D2Q9
    if fd==0
        % for loop
%         for l=1:q
%             f_eq(l,1)=w(l)*rho*(1+3*(V(:,l)')*U+9/2*((V(:,l)')*U)^2-3/2*(U')*U);
%         end
        % vectorized
        f_eq=rho*(1+3*V'*U+9/2*((V'*U).*(V'*U))-3/2*(U')*U).*w'; % Times
        %f_eq=diag(rho*(1+3*V'*U+9/2*diag((V'*U)*(V'*U)')-3/2*(U')*U)*w); % diag(matrix)
        
    elseif fd==1
        % for loop
%         for l=1:q
%             f_eq(l,1)=w(l)*(rho+Rho_r*(3*(V(:,l)')*U+9/2*((V(:,l)')*U)^2-3/2*(U')*U));
%         end
        % Vectorized
        f_eq=rho+Rho_r*(3*V'*U+9/2*((V'*U).*(V'*U))-3/2*(U')*U).*w'; % Times
        %f_eq=diag(rho+Rho_r*(3*V'*U+9/2*diag((V'*U)*(V'*U)')-3/2*(U')*U)*w); % diag(matrix)
    else
        error('Flag for density is invalid or unavailable!');
    end
elseif q==13 % D2Q13
    if fd==0
        % For loop
%         for l=1:q
%             f_eq(l,1)=w(l)*rho*(1+2*(V(:,l)')*U+2*((V(:,l)')*U)^2-(U')*U+4*((V(:,l)')*U)^3-6*((V(:,l)')*U)*((U')*U));
%         end
        % vectorized
        f_eq=diag(rho*(1+2*V'*U+2*diag((V'*U)*(V'*U)')-(U')*U+4*diag((V'*U)*(diag((V'*U)*(V'*U)'))')-6*(V'*U)*((U')*U))*w);
    elseif fd==1
        % For loop
%         for l=1:q
%             f_eq(l,1)=w(l)*(rho+Rho_r*(2*(V(:,l)')*U+2*((V(:,l)')*U)^2-(U')*U+4*((V(:,l)')*U)^3-6*((V(:,l)')*U)*((U')*U)));
%         end
        % vectorized
        f_eq=diag(rho+Rho_r*(2*V'*U+2*diag((V'*U)*(V'*U)')-(U')*U+4*diag((V'*U)*(diag((V'*U)*(V'*U)'))')-6*(V'*U)*((U')*U))*w);
    else
        error('Flag for density is invalid or unavailable!');
    end
elseif q==37 % D2Q37
    c_s=1; % fix speed of sound to be 1
    if fd==0
        f_eq=rho*w'.*(1+V'*U/c_s^2 ...
                       +((V'*U).*(V'*U)-c_s^2*(U')*U+c_s^2*(T-1)*(diag(V'*V)-2*c_s^2))/2/c_s^4 ...
                       +(V'*U).*((V'*U).*(V'*U)-3*c_s^2*((U')*U)+3*c_s^2*(T-1)*(diag(V'*V)-4*c_s^2))/6/c_s^6 ...
                       +(((V'*U).*(V'*U)).*((V'*U).*(V'*U))-6*c_s^2*((U')*U)*((V'*U).*(V'*U))+3*c_s^4*((U')*U)^2 ...
                       +6*c_s^2*(T-1)*(((V'*U).*(V'*U)).*(diag(V'*V)-6*c_s^2)+c_s^2*(U'*U)*(4*c_s^2-diag(V'*V))) ...
                       +3*c_s^4*(T-1)^2*(diag(V'*V).*diag(V'*V)-8*c_s^2*diag(V'*V)+8*c_s^4))/24/c_s^8); % Possibility 1
%         f_eq=rho*w'.*(1+V'*U/c_s^2 ...
%                        +((V'*U).*(V'*U)-c_s^2*(U')*U)/2/c_s^4 ...
%                        +(V'*U).*((V'*U).*(V'*U)-3*c_s^2*((U')*U))/6/c_s^6 ...
%                        +(((V'*U).*(V'*U)).*((V'*U).*(V'*U))-6*c_s^2*((V'*U).*(V'*U))+3*c_s^4*((U')*U)^2)/24/c_s^8); % Possibility 2, which is incorrect
                   
    elseif fd==1
        error('Flag for density is invalid or unavailable!');
    else
        error('Flag for density is invalid or unavailable!');
    end
else
    ;
end