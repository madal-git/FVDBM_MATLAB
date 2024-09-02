function bool=in_triangle_sub(P,T1,T2,T3)
% bool=in_triangle_sub(P,T1,T2,T3) returns 1 (yes) or no (0) that whether the
% point P is within the triangle with three vertices T1, T2 and T3
% All coordinate must be column vector
% This is a specific funtion used to assist on_edge and in_triangle

De=[1,T1(1,1),T1(2,1);1,T2(1,1),T2(2,1);1,T3(1,1),T3(2,1)];
Ta=[1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2];

De1=[1,T1(1,1),T1(2,1);1,T2(1,1),T2(2,1);1,P(1,1),P(2,1)];
Ta1=[1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2];

De2=[1,T2(1,1),T2(2,1);1,T3(1,1),T3(2,1);1,P(1,1),P(2,1)];
Ta2=[1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2;1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2];

De3=[1,T3(1,1),T3(2,1);1,T1(1,1),T1(2,1);1,P(1,1),P(2,1)];
Ta3=[1,T3(1,1),T3(2,1),T3(1,1)^2+T3(2,1)^2;1,T1(1,1),T1(2,1),T1(1,1)^2+T1(2,1)^2;1,P(1,1),P(2,1),P(1,1)^2+P(2,1)^2;1,T2(1,1),T2(2,1),T2(1,1)^2+T2(2,1)^2];

if det(De)*det(Ta)<0 && det(De1)*det(Ta1)>0 && det(De2)*det(Ta2)>0 && det(De3)*det(Ta3)>0
    bool=1;
else
    bool=0;
end