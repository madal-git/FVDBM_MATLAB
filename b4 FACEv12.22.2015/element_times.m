function u_times_v=element_times(u,v)

l1=length(u);
l2=length(v);
if l1~=l2
    error('The lengths of two vectors are different!');
else
    u_times_v=zeros(l1,1);
    for s=1:l1
        u_times_v(s)=u(s)*v(s);
    end
end