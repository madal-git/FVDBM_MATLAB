x_minmod_PC_PL=C_cut_x(1,:)
y_minmod_PC_PL=C_cut_y(2,:)
u_minmod_PC_PL=U_cut(1,:)
v_minmod_PC_PL=V_cut(2,:)
u_minmod_PC_PL=u_minmod_PC_PL/0.4
v_minmod_PC_PL=v_minmod_PC_PL/0.4
y_minmod_PC_PL=y_minmod_PC_PL/2
x_minmod_PC_PL=x_minmod_PC_PL/2

x_SUPERBEE_PC_PL=C_cut_x(1,:)
y_SUPERBEE_PC_PL=C_cut_y(2,:)
u_SUPERBEE_PC_PL=U_cut(1,:)
v_SUPERBEE_PC_PL=V_cut(2,:)
u_SUPERBEE_PC_PL=u_SUPERBEE_PC_PL/0.4
v_SUPERBEE_PC_PL=v_SUPERBEE_PC_PL/0.4
y_SUPERBEE_PC_PL=y_SUPERBEE_PC_PL/2
x_SUPERBEE_PC_PL=x_SUPERBEE_PC_PL/2


x_minmod_PC_PP=C_cut_x(1,:)
y_minmod_PC_PP=C_cut_y(2,:)
u_minmod_PC_PP=U_cut(1,:)
v_minmod_PC_PP=V_cut(2,:)
u_minmod_PC_PP=u_minmod_PC_PP/0.4
v_minmod_PC_PP=v_minmod_PC_PP/0.4
y_minmod_PC_PP=y_minmod_PC_PP/2
x_minmod_PC_PP=x_minmod_PC_PP/2

x_SUPERBEE_PC_PP=C_cut_x(1,:)
y_SUPERBEE_PC_PP=C_cut_y(2,:)
u_SUPERBEE_PC_PP=U_cut(1,:)
v_SUPERBEE_PC_PP=V_cut(2,:)
u_SUPERBEE_PC_PP=u_SUPERBEE_PC_PP/0.4
v_SUPERBEE_PC_PP=v_SUPERBEE_PC_PP/0.4
y_SUPERBEE_PC_PP=y_SUPERBEE_PC_PP/2
x_SUPERBEE_PC_PP=x_SUPERBEE_PC_PP/2


x_minmod_PL_PP=C_cut_x(1,:)
y_minmod_PL_PP=C_cut_y(2,:)
u_minmod_PL_PP=U_cut(1,:)
v_minmod_PL_PP=V_cut(2,:)
u_minmod_PL_PP=u_minmod_PL_PP/0.4
v_minmod_PL_PP=v_minmod_PL_PP/0.4
y_minmod_PL_PP=y_minmod_PL_PP/2
x_minmod_PL_PP=x_minmod_PL_PP/2

x_SUPERBEE_PL_PP=C_cut_x(1,:)
y_SUPERBEE_PL_PP=C_cut_y(2,:)
u_SUPERBEE_PL_PP=U_cut(1,:)
v_SUPERBEE_PL_PP=V_cut(2,:)
u_SUPERBEE_PL_PP=u_SUPERBEE_PL_PP/0.4
v_SUPERBEE_PL_PP=v_SUPERBEE_PL_PP/0.4
y_SUPERBEE_PL_PP=y_SUPERBEE_PL_PP/2
x_SUPERBEE_PL_PP=x_SUPERBEE_PL_PP/2

% x_PC=C_cut_x(1,:)
% y_PC=C_cut_y(2,:)
% u_PC=U_cut(1,:)
% v_PC=V_cut(2,:)
% u_PC=u_PC/0.4
% v_PC=v_PC/0.4
% y_PC=y_PC/2
% x_PC=x_PC/2
% 
% 
% x_PL=C_cut_x(1,:)
% y_PL=C_cut_y(2,:)
% u_PL=U_cut(1,:)
% v_PL=V_cut(2,:)
% u_PL=u_PL/0.4
% v_PL=v_PL/0.4
% y_PL=y_PL/2
% x_PL=x_PL/2
% 
% 
% x_PP=C_cut_x(1,:)
% y_PP=C_cut_y(2,:)
% u_PP=U_cut(1,:)
% v_PP=V_cut(2,:)
% u_PP=u_PP/0.4
% v_PP=v_PP/0.4
% y_PP=y_PP/2
% x_PP=x_PP/2
% 
% 
% 
% 
% 
u_sim=u_minmod_PC_PL;
y_sim=y_minmod_PC_PL;

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
        for j=1:H
            if y_Ghia(i)>y_sim(j) && y_Ghia(i)<y_sim(j+1)
                break;
            end
        end
        u_sim_sampled(i)=((y_sim_sampled(i)-y_sim(j))*u_sim(j+1)+(y_sim(j+1)-y_sim_sampled(i))*u_sim(j))/(y_sim(j+1)-y_sim(j));
    end
end

plot(u_sim,y_sim,u_sim_sampled,y_sim_sampled)
norm(u_Ghia-u_sim_sampled,2)


v_sim=v_SUPERBEE_PL_PP;
x_sim=x_SUPERBEE_PL_PP;

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
        for j=1:H
            if x_Ghia(i)>x_sim(j) && x_Ghia(i)<x_sim(j+1)
                break;
            end
        end
        v_sim_sampled(i)=((x_sim_sampled(i)-x_sim(j))*v_sim(j+1)+(x_sim(j+1)-x_sim_sampled(i))*v_sim(j))/(x_sim(j+1)-x_sim(j));
    end
end

plot(x_sim,v_sim,x_sim_sampled,v_sim_sampled)
norm(v_Ghia-v_sim_sampled,2)