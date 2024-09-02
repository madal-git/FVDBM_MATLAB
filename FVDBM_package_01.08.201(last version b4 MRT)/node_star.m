function x_node=node_star(ND,x,tag)
% x_node=node_star(ND,f) calculates the arbitrary vector at current node based
% on the vectors of the triangles connected to it
% ND is the NODE data structure of current node
% x is the matrix of triangle centroids of the entire domain
% tag is the flag that controls the scheme
% tag=0---normal scheme
% tag=1---periodic scheme

q=length(x(:,1));

if tag==0
    NP_NC=ND{5};
    NP_DC=ND{6};
    x_sum=zeros(q,1);
%     for m=1:q
%         for g=1:ND{4};
%             x_sum(m,1)=x_sum(m,1)+x(m,NP_NC(g))/NP_DC(g);
%         end
%     end;
    % Vectorzied
    for g=1:ND{4};
        x_sum(:,1)=x_sum(:,1)+x(:,NP_NC(g))/NP_DC(g);
    end
    
%     for m=1:q
%         x_node(m,1)=x_sum(m,1)/ND{7};
%     end
    % Vectorized
    x_node=x_sum(:,1)/ND{7};
elseif tag==1
    NP_NC=ND{9};
    NP_DC=ND{10};
    x_sum=zeros(q,1);
%     for m=1:q
%         for g=1:ND{12};
%             x_sum(m,1)=x_sum(m,1)+x(m,NP_NC(g))/NP_DC(g);
%         end
%     end;
    % Vectorized
    for g=1:ND{8};
        x_sum(:,1)=x_sum(:,1)+x(:,NP_NC(g))/NP_DC(g);
    end
    
%     for m=1:q
%         x_node(m,1)=x_sum(m,1)/ND{15};
%     end
    % Vectorized
    x_node=x_sum(:,1)/ND{11};
else
    error('The control for node_star() is unavailable or invalid!');
end