function bool=on_edge(P,ND1,ND2)
% function bool=on_edge(P,ND1,ND2) determines whether given point P is on
% the edges formed by linking ND1 and ND2, the other two points
% P is the coordinates of given point, has to be column vector
% ND1 is the coordinates of one end point of the edge, has to be column vector
% ND2 is the coordinates of one end point of the edge, has to be column vector
% bool is the boolen result. 1---The given point is on the edge; 0---No


% Algorithm 1: round-off error is 1e-2
% ref=10*dis(ND1,ND2);

% if double(ref+dis(P,ND1)+dis(P,ND2))==double(ref+dis(ND1,ND2))
%     bool=1;
% else
%     bool=0;
% end



% Algorithm 2: round-off error is 1e-8
% Fac=10;
% if double(Fac*(dis(P,ND1)+dis(P,ND2)))==double(Fac*dis(ND1,ND2))
%     bool=1;
% else
%     bool=0;
% end



% Algorithm 3: round-off error is 1e-50
if in_triangle_sub(single((P+ND1+ND2)/3),single(P),single(ND1),single(ND2))
% if in_triangle_sub((P+ND1+ND2)/3,P,ND1,ND2)
% if in_triangle_sub(double((P+ND1+ND2)/3),double(P),double(ND1),double(ND2))
    bool=0;
else
    bool=1;
end