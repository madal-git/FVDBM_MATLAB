N=1000;
tic;
for i=1:N
    f=f_old.*f_old;
    for j=1:100
    f_face=f_old+(f_old-abs(f_old).*f_old*0.32/2).*(f_old+f_old);
    end
end
time_per_iteration=toc/N


tic;
for i=1:N
    f=diag(f_old*f_old');
    for j=1:100
    f_face=f_old+diag((f_old-diag(abs(f_old)*f_old')*0.32/2)*(f_old+f_old)');
    end
end
time_per_iteration=toc/N

tic;
for i=1:N
    f=element_times(f_old,f_old);
    for j=1:100
    f_face=f_old+element_times(f_old-element_times(abs(f_old),f_old)*0.32/2,f_old+f_old);
    end
end
time_per_iteration=toc/N