dt_over_dx=0.05:0.05:0.5;
dt_over_tau=0.5:0.5:5;
L=length(dt_over_dx);
H=length(dt_over_tau);
A=cell(H,L);
%% Initialization and Default value for instability 
for i=1:L
    for j=1:H
        A{j,i}=[dt_over_tau(j), dt_over_dx(i), 0];
    end
end

%% Set stability value
for i=1:L
    for j=1:H
        s=A{j,i};
        if single(s(2))==single(0.05)
            if s(1)<2.5
                s(3)=1;
                A{j,i}=s;
            end
        end
        if single(s(2))==single(0.1)
            if s(1)<2
                s(3)=1;
                A{j,i}=s;
            end
        end
        if single(s(2))==single(0.15)
            if s(1)<1
                s(3)=1;
                A{j,i}=s;
            end
        end
%         if single(s(2))==single(0.2)
%             if s(1)<1.5
%                 s(3)=1;
%                 A{j,i}=s;
%             end
%         end
%         if single(s(2))==single(0.25)
%             if s(1)<1.5
%                 s(3)=1;
%                 A{j,i}=s;
%             end
%         end
%         if single(s(2))==single(0.3)
%             if s(1)<1.5
%                 s(3)=1;
%                 A{j,i}=s;
%             end
%         end
%         if single(s(2))==single(0.35)
%             if s(1)<1
%                 s(3)=1;
%                 A{j,i}=s;
%             end
%         end
    end
end


figure;
for i=1:L
    for j=1:H
        s=A{j,i};
        if s(3)==0
            plot(s(2),s(1),'Marker', 'x','Markersize',5, 'MarkerFaceColor', 'r','MarkerEdgeColor', 'r');;
        else
            plot(s(2),s(1),'Marker', 'o','Markersize',5, 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k');
        end
        hold on
    end
end
xlim([0 0.55])
hold on
ylim([0 5.5])
grid on
title('(d) RK4');
ax = gca;
ax.XTick = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55];
xlabel('\Deltat/\Deltax');
ylabel('\Deltat/\lambda');