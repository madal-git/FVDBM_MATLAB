u_sim=u;
y_sim=y;

L=length(u_Ghia);
H=length(y_sim);
u_sim_sampled=zeros(L,1);
y_sim_sampled=zeros(L,1);
for i=1:L
    y_sim_sampled(i)=y_Ghia(i);
    if i==1
        u_sim_sampled(i)=u_sim(end);
    elseif i==L
        u_sim_sampled(i)=u_sim(1);
    else
        for j=1:H-1
            if y_Ghia(i)>=y_sim(j) && y_Ghia(i)<=y_sim(j+1)
                break;
            end
        end
        u_sim_sampled(i)=((y_sim_sampled(i)-y_sim(j))*u_sim(j+1)+(y_sim(j+1)-y_sim_sampled(i))*u_sim(j))/(y_sim(j+1)-y_sim(j));
    end
end

figure
plot(u_sim,y_sim,u_sim_sampled,y_sim_sampled,u_Ghia,y_Ghia)
u_err_L2=norm(u_Ghia-u_sim_sampled,2)/norm(u_Ghia,2)


v_sim=v;
x_sim=x;

L=length(v_Ghia);
H=length(x_sim);
v_sim_sampled=zeros(L,1);
x_sim_sampled=zeros(L,1);
for i=1:L
    x_sim_sampled(i)=x_Ghia(i);
    if i==1
        v_sim_sampled(i)=v_sim(end);
    elseif i==L
        v_sim_sampled(i)=v_sim(1);
    else
        for j=1:H-1
            if x_Ghia(i)>=x_sim(j) && x_Ghia(i)<=x_sim(j+1)
                break;
            end
        end
        v_sim_sampled(i)=((x_sim_sampled(i)-x_sim(j))*v_sim(j+1)+(x_sim(j+1)-x_sim_sampled(i))*v_sim(j))/(x_sim(j+1)-x_sim(j));
    end
end

figure
plot(x_sim,v_sim,x_sim_sampled,v_sim_sampled,x_Ghia,v_Ghia)
v_err_L2=norm(v_Ghia-v_sim_sampled,2)/norm(v_Ghia,2)