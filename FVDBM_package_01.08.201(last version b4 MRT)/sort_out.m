function A=sort_out(a)
% function A=sort_out(a) picks out a series consisting of all unique values
% in a row vector a
% Check whether a is a row vector
s=size(a);
if s(1)~=1
    error('Input must be a row vector');
end
M=length(a);
for i=1:M
    a_temp=a(i+1:end);
    N=length(a_temp);
    counter=0;
    b=Inf;
    for j=1:N
        if a(i)~=a_temp(j)
            counter=counter+1;
            b(counter)=a_temp(j);
        end
    end
    if counter==0 % b==Inf
        a=a(1:i);
        break;
    else
        a=[a(1:i),b];
    end
    M=length(a);
end
A=a;