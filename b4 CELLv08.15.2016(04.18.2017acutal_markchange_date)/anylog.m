function A=anylog(B,C)
% A=anylog(B,C) returns the power A that satisfied
% B^A=C

if B<=0
    error('The base number must be positive');
end

if C<=0
    error('The result must be positive');
end

A_low=-1;
A_up=1;
while B^A_low>C
    A_low=A_low*10;
end
while B^A_up<C
    A_up=A_up*10;
end

A_temp=(A_low+A_up)/2;

while abs(C-B^A_temp)>1e-8
    if B^A_temp<C
        A_low=A_temp;
    else
        A_up=A_temp;
    end
    A_temp=(A_low+A_up)/2;
end
A=A_temp;