T_max=1;
L=2;
N=128;
Alpha=0.002;
dL=L/(N-1);
x=0:dL:L;
y=0:dL:L;
T=zeros(N,N);
T_new=zeros(N,N);
T_ana=zeros(N,N);
Timer=0;
Time_stop=48.828125; %% seconds
dt_max=(dL*(N-1)/N)^2/Alpha/2;
dt=dt_max/2;
%% Initial Condition
for i=1:N
    for j=1:N
        T(i,j)=T_max*exp(-((x(i)-L/2)^2+(y(j)-L/2)^2)/(L/2)^2);
% if i>=N/4 && i<=3*N/4 && j>=N/4 && j<=3*N/4
%     T(i,j)=T_max;
% end
    end
end
% Finite Difference solution
while Timer<Time_stop
    for i=1:N
        for j=1:N
            if i==1
                if j==1
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N-1,j)+T(i,j+1)+T(i,N-1));
                elseif j==N
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N-1,j)+T(i,2)+T(i,j-1));
                else
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N-1,j)+T(i,j+1)+T(i,j-1));
                end
            elseif i==N
                if j==1
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(2,j)+T(i-1,j)+T(i,j+1)+T(i,N-1));
                elseif j==N
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(2,j)+T(i-1,j)+T(i,2)+T(i,j-1));
                else
                    T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(2,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
                end
            elseif j==1
                T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,N-1));
            elseif j==N
                T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,2)+T(i,j-1));
            else
                T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
            end
        end
    end
    T=T_new;
%     for i=1:N
%         for j=1:N
%             if i==1
%                 if j==1
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N,j)+T(i,j+1)+T(i,N));
%                 elseif j==N
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N,j)+T(i,1)+T(i,j-1));
%                 else
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(N,j)+T(i,j+1)+T(i,j-1));
%                 end
%             elseif i==N
%                 if j==1
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(1,j)+T(i-1,j)+T(i,j+1)+T(i,N));
%                 elseif j==N
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(1,j)+T(i-1,j)+T(i,1)+T(i,j-1));
%                 else
%                     T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
%                 end
%             elseif j==1
%                 T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,N));
%             elseif j==N
%                 T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,1)+T(i,j-1));
%             else
%                 T_new(i,j)=(1-4*dt*Alpha/dL^2)*T(i,j)+dt*Alpha/dL^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
%             end
%         end
%     end
%     T=T_new;
    Timer=Timer+dt;
end
figure
contourf(T_new,100);
axis equal tight;
%% Analytical solution
for i=1:N
    for j=1:N
        T_ana(i,j)=T_max/sqrt(1+4*(Timer-1)*dt*Alpha/(L/2)^2)*exp(-((x(i)-L/2)^2+(y(j)-L/2)^2)/((L/2)^2+4*(Timer-1)*dt*Alpha));
    end
end
% hold on
% figure
% surf(T_ana)

%%Sampling the data at specified locations in order to compare with results
%%from FVDBM
M=37; % The sampling points are discrete locations at the crosspoints of M-by-M grid
T_ana_sampled=zeros(M,M);
X=0:L/(M-1):L;
Y=0:L/(M-1):L;
for p=1:M
    for q=1:M
        if p==1
            if q==1
                T_ana_sampled(p,q)=T_new(1,1);
            elseif q==M
                T_ana_sampled(p,q)=T_new(1,N);
            else
                for s=1:N-1
                    if Y(q)>=y(s) && Y(q)<=y(s+1)
                        T_ana_sampled(p,q)=(Y(q)-y(s))/dL*T_new(1,s+1)+(y(s+1)-Y(q))/dL*T_new(1,s);
                    end
                end
            end
        elseif p==M
            if q==1
                T_ana_sampled(p,q)=T_new(N,1);
            elseif q==M
                T_ana_sampled(p,q)=T_new(N,N);
            else
                for s=1:N-1
                    if Y(q)>=y(s) && Y(q)<=y(s+1)
                        T_ana_sampled(p,q)=(Y(q)-y(s))/dL*T_new(N,s+1)+(y(s+1)-Y(q))/dL*T_new(N,s);
                    end
                end
            end
        elseif q==1
            for s=1:N-1
                if X(p)>=x(s) && X(p)<=x(s+1)
                    T_ana_sampled(p,q)=(X(p)-x(s))/dL*T_new(s+1,1)+(x(s+1)-X(p))/dL*T_new(s,1);
                end
            end
        elseif q==M
            for s=1:N-1
                if X(p)>=x(s) && X(p)<=x(s+1)
                    T_ana_sampled(p,q)=(X(p)-x(s))/dL*T_new(s+1,N)+(x(s+1)-X(p))/dL*T_new(s,N);
                end
            end
        else
            for s=1:N-1
                for t=1:N-1
                    if (X(p)>=x(s) && X(p)<=x(s+1)) && (Y(q)>=y(t) && Y(q)<=y(t+1))
                        if single(L+X(p))==single(L+x(s)) && single(L+Y(q))==single(L+y(t))
                            T_ana_sampled(p,q)=T_new(s,t);
                        elseif single(L+X(p))==single(L+x(s)) && single(L+Y(q))==single(L+y(t+1))
                            T_ana_sampled(p,q)=T_new(s,t+1);
                        elseif single(L+X(p))==single(L+x(s+1)) && single(L+Y(q))==single(L+y(t))
                            T_ana_sampled(p,q)=T_new(s+1,t);
                        elseif single(L+X(p))==single(L+x(s+1)) && single(L+Y(q))==single(L+y(t+1))
                            T_ana_sampled(p,q)=T_new(s+1,t+1);
                        else
                            Quad_1=[x(s);y(t)];
                            Quad_2=[x(s);y(t+1)];
                            Quad_3=[x(s+1);y(t)];
                            Quad_4=[x(s+1);y(t+1)];
                            Dis_1=dis(Quad_1,[X(p);Y(q)]);
                            Dis_2=dis(Quad_2,[X(p);Y(q)]);
                            Dis_3=dis(Quad_3,[X(p);Y(q)]);
                            Dis_4=dis(Quad_4,[X(p);Y(q)]);
                            T_ana_sampled(p,q)=(T_new(s,t)/Dis_1+T_new(s,t+1)/Dis_2+T_new(s+1,t)/Dis_3+T_new(s+1,t+1)/Dis_4)/(1/Dis_1+1/Dis_2+1/Dis_3+1/Dis_4);
                        end
                    end
                end
            end
        end
    end
end
figure;
contourf(T_ana_sampled,100);
axis equal tight;
figure
surf(ones(N,N).*x,(ones(N,N).*(y))',T_new);
hold on
surf(ones(M,M).*X,(ones(M,M).*Y)',T_ana_sampled);
sum(sum(T_new))/N^2
sum(sum(T_ana_sampled))/M^2
if M==19
    save('T_ana_IRT18by18.mat','T_ana_sampled');
elseif M==37
    save('T_ana_IRT36by36.mat','T_ana_sampled');
elseif M==73
    save('T_ana_IRT72by72.mat','T_ana_sampled');
else
    error('Incorrect dimension!');
end