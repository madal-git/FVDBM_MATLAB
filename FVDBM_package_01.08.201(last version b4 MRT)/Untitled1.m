n=length(V(1,:));
Ori=[0;0];
for i=1:n
    if n==37
    plot([Ori(1),V(1,i)]/c_l,[[Ori(2),V(2,i)]]/c_l);
    else
        plot([Ori(1),V(1,i)],[[Ori(2),V(2,i)]]);
    end
    hold on
end
axis equal
grid on