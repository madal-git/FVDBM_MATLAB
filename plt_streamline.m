if strcmp(top,'Moving Wall')==1 && strcmp(right,'Stationary Wall')==1 && strcmp(bottom,'Stationary Wall')==1 && strcmp(left,'Stationary Wall')==1 % Lid-driven cavity flow
    streamline_ldcf(U_m(1,1),N_L,N_H,N_I,X1,X2,Y1,Y2,U,U_nd,CELL,M,NODE,N)
else
    error('The streamline plot cannot be generated for the current flow type!');
end