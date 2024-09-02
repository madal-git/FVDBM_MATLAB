function fl=flux_check(P,f_nd1,f_nd2,f_nd3,V,e,q)
% fl=flux_check(P,f_nd,V,e,q) calculates flux 
% for single edge and specific velocity direction of 
% individual triangle.
% P is the current triangle
% f_nd1 is the pdf vector at node 1 of current triangle
% f_nd2 is the pdf vector at node 2 of current triangle
% f_nd3 is the pdf vector at node 3 of current triangle
% V is the lattice velocity matrix
% e is numbering of edge of current triangle
%%%% 1----first edge; 2----second edge; 3----third edge
% g is the numbering of lattice velocity component corresponding to V
% fl is the scalar flux of pdf direction q though the edge e of current triangle


switch e
    case 1
        fl=P{34}*V(:,q)*(f_nd1(q,1)+f_nd2(q,1))/2*P{31};
    case 2
        fl=P{35}*V(:,q)*(f_nd2(q,1)+f_nd3(q,1))/2*P{32};
    case 3
        fl=P{36}*V(:,q)*(f_nd3(q,1)+f_nd1(q,1))/2*P{33};
end