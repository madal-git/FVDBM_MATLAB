function v_o=norm_orthog(P1,P2)
% v_o=norm_orthog(P1,P2) returns a unit vertor that is orthogonal to the
% vector that connects two points, P1 and P2

% P1 and P2 have to be column vectors


[D1,L1]=size(P1);
[D2,L2]=size(P1);

if D1~=D2 || L1~=L2
    error('Check the dimensions of the coordinates of input points');
end

if D1~=2
    error('The current function can only deal with 2D!');
end

if L1~=1
    error('The input is not a point');
end

v_o=norm_v([-(P2(2,1)-P1(2,1));(P2(1,1)-P1(1,1))]);