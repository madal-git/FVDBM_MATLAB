function bool=in_triangle(P,T1,T2,T3)
% bool=in_triangle(P,T1,T2,T3) returns 1 (yes) or no (0) that whether the
% point P is within the triangle with three vertices T1, T2 and T3
% All coordinate must be column vector

%% ALGORITHM 1: Sequence of T1, T2 and T3 matter
% De=[1,T1(1,1),T1(2,1);1,T2(1,1),T2(2,1);1,T3(1,1),T3(2,1)];
% Ta=[1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2];
% 
% De1=[1,T1(1,1),T1(2,1);1,T2(1,1),T2(2,1);1,P(1,1),P(2,1)];
% Ta1=[1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2];
% 
% De2=[1,T2(1,1),T2(2,1);1,T3(1,1),T3(2,1);1,P(1,1),P(2,1)];
% Ta2=[1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2];
% 
% De3=[1,T3(1,1),T3(2,1);1,T1(1,1),T1(2,1);1,P(1,1),P(2,1)];
% Ta3=[1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2];
% 
% if det(De)*det(Ta)<0 && det(De1)*det(Ta1)>0 && det(De2)*det(Ta2)>0 && det(De3)*det(Ta3)>0
%     if on_edge(P,T1,T2) || (on_edge(P,T2,T3) || on_edge(P,T3,T1))
%         bool=0;
%     else
%         bool=1;
%     end
% else
%     bool=0;
% end

% %% ALGORITHM 2: Sequence of T1, T2 and T3 matter does not matter
% if in_triangle_sub(P,T1,T2,T3)
%     if (on_edge(P,T1,T2) || (on_edge(P,T2,T3) || on_edge(P,T3,T1))) || (on_edge(P,T2,T1) || (on_edge(P,T3,T2) || on_edge(P,T1,T3)))
%         bool=0;
%     else
%         bool=1;
%     end
% else
%     bool=0;
% end

%% ALGORITHM 3: Dilation, shrink the ones over higher limit, enlarge the ones under lower limits
% L_th_low=0.1;
% L_th_high=0.2;
% L=(dis(T1,T2)+dis(T2,T3)+dis(T3,T1))/3;
% if L>=L_th_low && L<=L_th_high
%     if in_triangle_sub(P,T1,T2,T3)
%         if on_edge(P,T1,T2) || (on_edge(P,T2,T3) || on_edge(P,T3,T1))
%             bool=0;
%         else
%             bool=1;
%         end
%     else
%         bool=0;
%     end
% else
%     if L<L_th_low
%         D_factor=10^ceil(log10(L_th_low/L));
%     else
%         D_factor=10^(-ceil(log10(L/L_th_high)));
%     end
%     Ori=(T1+T2+T3)/3;
%     P_new=(P-Ori)*D_factor;
%     T1_new=(T1-Ori)*D_factor;
%     T2_new=(T2-Ori)*D_factor;
%     T3_new=(T3-Ori)*D_factor;
%     if in_triangle_sub(P_new,T1_new,T2_new,T3_new)
%         if on_edge(P_new,T1_new,T2_new) || (on_edge(P_new,T2_new,T3_new) || on_edge(P_new,T3_new,T1_new))
%             bool=0;
%         else
%             bool=1;
%         end
%     else
%         bool=0;
%     end
% end


%% ALGORITHM 4: Dilation, Change all triangle size to target size
L_th_target=0.04;
L=(dis(T1,T2)+dis(T2,T3)+dis(T3,T1))/3;
D_factor=L_th_target/L;
Ori=(T1+T2+T3)/3;
P_new=(P-Ori)*D_factor;
T1_new=(T1-Ori)*D_factor;
T2_new=(T2-Ori)*D_factor;
T3_new=(T3-Ori)*D_factor;
if in_triangle_sub(P_new,T1_new,T2_new,T3_new)
    if on_edge(P_new,T1_new,T2_new) || (on_edge(P_new,T2_new,T3_new) || on_edge(P_new,T3_new,T1_new))
        bool=0;
    else
        bool=1;
    end
else
    bool=0;
end