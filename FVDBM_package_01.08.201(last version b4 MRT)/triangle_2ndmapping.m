function x=triangle_2ndmapping(C1,C2,C3,x1,x2,x3,C)
% function x=triangle_2ndmapping(C1,C2,C3,x1,x2,x3,C) calculates the value
% of x, which is located at point C in a linear plane bound by point C1,
% C2 and C3 with unknown values x1, x2 and x3 at each point.
% C1 is the coordinate of first point, column vector
% C2 is the coordinate of second point, column vector
% C3 is the coordinate of third point, column vector
% x1 is the value at point C1, a scalar or row vector
% x2 is the value at point C2, a scalar or row vector
% x3 is the value at point C3, a scalar or row vector
% C is the coordinate of caculated point, column vector
% x is the linear interpolated value at point C,a row vector

A=[C1(1,1),C1(2,1),1;C2(1,1),C2(2,1),1;C3(1,1),C3(2,1),1];
F=[x1;x2;x3];
B=A\F;
x=[C(1,1),C(2,1),1]*B;
% %% 2nd method, if x1, x2 and x3 are column vectors
% A=[C1(1,1),C1(2,1),1;C2(1,1),C2(2,1),1;C3(1,1),C3(2,1),1];
% A_inverse=A\eye(size(A));
% F=[x1,x2,x3];
% B=F*A_inverse;
% x=F*[C(1,1);C(2,1);1];

% %% Third method
% x=[C(1,1),C(2,1),1]*([C1(1,1),C1(2,1),1;C2(1,1),C2(2,1),1;C3(1,1),C3(2,1),1]\[x1;x2;x3]);


%% Check
% x1_new=[C1(1,1),C1(2,1),1]*B;
% x2_new=[C2(1,1),C2(2,1),1]*B;
% x3_new=[C3(1,1),C3(2,1),1]*B;
% e=100;
% if single(e+norm(x1-x1_new,1))~=single(e) || single(e+norm(x2-x2_new,1))~=single(e) || single(e+norm(x3-x3_new,1))~=single(e)
%     error('wrong!');
% end