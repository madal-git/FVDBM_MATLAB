function [matrix,matrix_inverse] = mrt_2d_trans_matrix (d,q)
%% function [matrix,matrix_inverse] = mrt_2d_trans_matrix (d,q) generetes the transformation matrix and its inverse matrix for the MRT model
%  d is the number of dimensions
%  q is the number of velocity components in the velocity space
%  matrix is the transformation matrix
%  matrix_inverse is the inverse of the transformation matrix

if d==1
    error('Temporarily not available!');
elseif d==2
    if q==5
        error('Temporarily not available!');
    elseif q==7
        error('Temporarily not available!');
    elseif q==9
        matrix=[1  1  1  1  1  1  1  1  1;...
               -4 -1 -1 -1 -1  2  2  2  2;...
                4 -2 -2 -2 -2  1  1  1  1;...
                0  1  0 -1  0  1 -1 -1  1;...
                0 -2  0  2  0  1 -1 -1  1;...
                0  0  1  0 -1  1  1 -1 -1;...
                0  0 -2  0  2  1  1 -1 -1;...
                0  1 -1  1 -1  0  0  0  0;...
                0  0  0  0  0  1 -1  1 -1];
         matrix_inverse=inv(matrix);
    else
        error('Wrong lattice model');
    end
elseif d==3
    error('Temporarily not available!');
else
    error('Wrong dimensions!');
end