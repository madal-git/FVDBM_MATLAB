function fcol=BGK(q,P,w,tau,f,f_nd1,f_nd2,f_nd3,f_eq,f_eq_nd1,f_eq_nd2,f_eq_nd3)
% fcol=BGK(q,P,w,tau,f,f_nd1,f_nd2,f_nd3,f_eq,f_eq_nd1,f_eq_nd2,f_eq_nd3)
% calculates the integral post-collision pdf on current 
% triangle by BGK collision model
% q is the total number of velocity components in the lattice
% P is the current triangle
% w is the weighting factor
% tau is the relaxation time
% f is the pdf vectior at the centroid of current triangle
% f_nd1 is the pdf vectior at the first node of current triangle
% f_nd2 is the pdf vectior at the second node of current triangle
% f_nd3 is the pdf vectior at the third node of current triangle
% f_eq is the equilibrium pdf vectior at the centroid of current triangle
% f_eq_nd1 is the equilibrium pdf vectior at the first node of current triangle
% f_eq_nd2 is the equilibrium pdf vectior at the second node of current triangle
% f_eq_nd3 is the equilibrium pdf vectior at the third node of current triangle
% fcol is the after-collision pdf vectior at the centroid of current triangle


fcol(:,1)=(-1)*P{6}/tau*(f(:,1)-f_eq(:,1));
%fcol(:,1)=(-1)*P{6}/tau*(w*(f(:,1)-f_eq(:,1))+(1-w)/3*((f_nd1(:,1)-f_eq_nd1(:,1))+(f_nd2(:,1)-f_eq_nd2(:,1))+(f_nd3(:,1)-f_eq_nd3(:,1))));