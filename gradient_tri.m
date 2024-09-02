function grd = gradient_tri(q,M,CELL,x_cell,x_node,fpdc)
%% function grd = gradient_tri(q,M,CELL,x_cell,x_node) calculates the spacial gradient of variable x
%  at the centroids of triangular mesh with the known values of x at the
%  centroids (x_cell) and the known values of x at the nodes (x_node)

% q is the dimension of x_cell and x_node
% M is the number of triangular cells
% CELL is the CELL data structure
% x_cell is the values of x at the centroids of all of the triangular cells
% x_node is the values of x at the nodes of the triangular mesh
% fpdc is the flag for the periodic conditions of the boundaries
grd=zeros(2,q,M);

%% Algorithm 1 see paper Maik Stiebler, Jonas Tolke, Manfred Krafczyk. An upwind discretization scheme for the finite volume lattice. Computers & Fluids
for k=1:M
    CL=CELL{k};
    B=zeros(2,q);
    A=0;
    if CL{37}==1 % Cell type 1
        cell_neigh=CL{38};
        w=CL{39};
        dis_vec=[CL{40},CL{41},CL{42}];
        A=CL{43};
        for i=1:q
            for j=1:3
                B(:,i)=B(:,i)+w(j)*dis_vec(:,j)*(x_cell(i,cell_neigh(j))-x_cell(i,k));
            end
        end
    elseif CL{37}==3
        cell_neigh=CL{44};
        nd=CL{45};
        w=CL{46};
        dis_vec=[CL{47},CL{48},CL{49}];
        A=CL{50};
        if cell_neigh(1)==0
            f_neigh=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
        elseif cell_neigh(2)==0
            f_neigh=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
        elseif cell_neigh(3)==0
            f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
        else
            error('logic error!');
        end
        for i=1:q
            for j=1:3
                B(:,i)=B(:,i)+w(j)*dis_vec(:,j)*(f_neigh(i,j)-x_cell(i,k));
            end
        end
    elseif CL{37}==2
        if fpdc==0 % No periodic boundaries
            cell_neigh=CL{44};
            nd=CL{45};
            w=CL{46};
            dis_vec=[CL{47},CL{48},CL{49}];
            A=CL{50};
        elseif fpdc==1 % Left and right boundaries are periodic
            cell_neigh=CL{58};
            nd=CL{59};
            w=CL{60};
            dis_vec=[CL{61},CL{62},CL{63}];
            A=CL{64};
        elseif fpdc==2 % Top and bottom boundaries are periodic
            cell_neigh=CL{65};
            nd=CL{66};
            w=CL{67};
            dis_vec=[CL{68},CL{69},CL{70}];
            A=CL{71};
        elseif fpdc==3 % All outer boundaries are periodic
            cell_neigh=CL{51};
            nd=CL{52};
            w=CL{53};
            dis_vec=[CL{54},CL{55},CL{56}];
            A=CL{57};
        else
            error('Wrong flag for periodic boundary conditions!');
        end
        len=length(setxor(cell_neigh,0));
        if len==4 % all neighbors are cells
            f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
        elseif len==2 % one neighbor cell is missing
            if cell_neigh(1)==0
                f_neigh=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
            elseif cell_neigh(2)==0
                f_neigh=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
            elseif cell_neigh(3)==0
                f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
            else
                error('logic error!');
            end
        else
            error('logic error!');
        end
        for i=1:q
            for j=1:3
                B(:,i)=B(:,i)+w(j)*dis_vec(:,j)*(f_neigh(i,j)-x_cell(i,k));
            end
        end
    else
        error('Logic error!');
    end
    
    for i=1:q
        grd(:,i,k)=A\B(:,i);
    end
end

  %% Algorithm 2 see paper W. Kyle Anderson and Daryl L. Bonhaus. An Implicit Upwind Algorithm for Computing Turbulent Flows on Unstructured Grids. Computers & Fluids
% for k=1:M
%     CL=CELL{k};
%     g_cent_x=zeros(q,1);
%     g_cent_y=zeros(q,1);
%     if CL{37}==1 % Cell type 1
%         cell_neigh=CL{38};
%         w_x=CL{39};
%         w_y=CL{40};
%         x_cent=x_cell(:,k);
%         x_sate=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%         for j=1:3
%             g_cent_x=g_cent_x+w_x(1,j)*(x_sate(:,j)-x_cent);
%             g_cent_y=g_cent_y+w_y(1,j)*(x_sate(:,j)-x_cent);
%         end
%     elseif CL{37}==3
%         cell_neigh=CL{44};
%         nd=CL{45};
%         w_x=CL{46};
%         w_y=CL{47};
%         x_cent=x_cell(:,k);
%         if cell_neigh(1)==0
%             x_sate=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%         elseif cell_neigh(2)==0
%             x_sate=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
%         elseif cell_neigh(3)==0
%             x_sate=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
%         else
%             error('logic error!');
%         end
%         for j=1:3
%             g_cent_x=g_cent_x+w_x(1,j)*(x_sate(:,j)-x_cent);
%             g_cent_y=g_cent_y+w_y(1,j)*(x_sate(:,j)-x_cent);
%         end
%     elseif CL{37}==2
%         if fpdc==0 % No periodic boundaries
%             cell_neigh=CL{44};
%             nd=CL{45};
%             w_x=CL{46};
%             w_y=CL{47};
%         elseif fpdc==1 % Left and right boundaries are periodic
%             cell_neigh=CL{58};
%             nd=CL{59};
%             w_x=CL{60};
%             w_y=CL{61};
%         elseif fpdc==2 % Top and bottom boundaries are periodic
%             cell_neigh=CL{65};
%             nd=CL{66};
%             w_x=CL{67};
%             w_y=CL{68};
%         elseif fpdc==3 % All outer boundaries are periodic
%             cell_neigh=CL{51};
%             nd=CL{52};
%             w_x=CL{53};
%             w_y=CL{54};
%         else
%             error('Wrong flag for periodic boundary conditions!');
%         end
%         len=length(setxor(cell_neigh,0));
%         if len==4 % all neighbors are cells
%             x_sate=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%         elseif len==2 % one neighbor cell is missing
%             if cell_neigh(1)==0
%                 x_sate=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%             elseif cell_neigh(2)==0
%                 x_sate=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
%             elseif cell_neigh(3)==0
%                 x_sate=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
%             else
%                 error('logic error!');
%             end
%         else
%             error('logic error!');
%         end
%         for j=1:3
%             g_cent_x=g_cent_x+w_x(1,j)*(x_sate(:,j)-x_cent);
%             g_cent_y=g_cent_y+w_y(1,j)*(x_sate(:,j)-x_cent);
%         end
%     else
%         error('Logic error!');
%     end
%     
%     grd(1,:,k)=g_cent_x';
%     grd(2,:,k)=g_cent_y';
% end

%% Algorithm 3 weighted scheme based on the distance to the center centroid
% for k=1:M
%     CL=CELL{k};
%     if CL{37}==1 % Cell type 1
%         cell_neigh=CL{38};
%         w=CL{39};
%         dis_vec=[CL{40},CL{41},CL{42}];
%         A=zeros(2,q);
%         for j=1:3
%             grd(:,:,k)=A+w(j)*((1./dis_vec(:,j))*(x_cell(:,cell_neigh(j))-x_cell(:,k))');
%         end
%     elseif CL{37}==3
%         cell_neigh=CL{44};
%         nd=CL{45};
%         w=CL{46};
%         dis_vec=[CL{47},CL{48},CL{49}];
%         A=zeros(2,q);
%         if cell_neigh(1)==0
%             f_neigh=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%         elseif cell_neigh(2)==0
%             f_neigh=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
%         elseif cell_neigh(3)==0
%             f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
%         else
%             error('logic error!');
%         end
%         for j=1:3
%             grd(:,:,k)=A+w(j)*((1./dis_vec(:,j))*(f_neigh(:,j)-x_cell(:,k))');
%         end
%     elseif CL{37}==2
%         if fpdc==0 % No periodic boundaries
%             cell_neigh=CL{44};
%             nd=CL{45};
%             w=CL{46};
%             dis_vec=[CL{47},CL{48},CL{49}];
%         elseif fpdc==1 % Left and right boundaries are periodic
%             cell_neigh=CL{58};
%             nd=CL{59};
%             w=CL{60};
%             dis_vec=[CL{61},CL{62},CL{63}];
%         elseif fpdc==2 % Top and bottom boundaries are periodic
%             cell_neigh=CL{65};
%             nd=CL{66};
%             w=CL{67};
%             dis_vec=[CL{68},CL{69},CL{70}];
%         elseif fpdc==3 % All outer boundaries are periodic
%             cell_neigh=CL{51};
%             nd=CL{52};
%             w=CL{53};
%             dis_vec=[CL{54},CL{55},CL{56}];
%         else
%             error('Wrong flag for periodic boundary conditions!');
%         end
%         len=length(setxor(cell_neigh,0));
%         if len==4 % all neighbors are cells
%             f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%         elseif len==2 % one neighbor cell is missing
%             if cell_neigh(1)==0
%                 f_neigh=[(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(2)),x_cell(:,cell_neigh(3))];
%             elseif cell_neigh(2)==0
%                 f_neigh=[x_cell(:,cell_neigh(1)),(x_node(:,nd(1))+x_node(:,nd(2)))/2,x_cell(:,cell_neigh(3))];
%             elseif cell_neigh(3)==0
%                 f_neigh=[x_cell(:,cell_neigh(1)),x_cell(:,cell_neigh(2)),(x_node(:,nd(1))+x_node(:,nd(2)))/2];
%             else
%                 error('logic error!');
%             end
%         else
%             error('logic error!');
%         end
%         A=zeros(2,q);
%         for j=1:3
%             grd(:,:,k)=A+w(j)*((1./dis_vec(:,j))*(f_neigh(:,j)-x_cell(:,k))');
%         end
%     else
%         error('Logic error!');
%     end
% end