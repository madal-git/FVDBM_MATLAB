    %% Filling FC{16}, {17}
    for l=1:O
        FC=FACE{l};
        S_V1=cell(K,1); % stencil for V1
        S_V2=cell(K,1); % stencil for V2
        neigh_up=FC{12};
        neigh_down=FC{13};
        %% S{1} and S{2}
        % S{1}, S{2} for V1
        DC1=zeros(1,q1);
        UC1=zeros(1,q1);
        for k=1:q1
            if FC{4}*V1(:,k)>=0
                DC1(1,k)=neigh_down(1,1);
                UC1(1,k)=neigh_up(1,1);
            else
                DC1(1,k)=neigh_up(1,1);
                UC1(1,k)=neigh_down(1,1);
            end
        end
        % Check
        % Cell pair
        for k=1:q1
            if length(unique(union([DC1(1,k),UC1(1,k)],[neigh_down(1,1),neigh_up(1,1)])))~=length([neigh_down(1,1),neigh_up(1,1)])
                error('The downwind & upwind cells are not correct!')
            end
        end
        % Direction
        for k=1:q1
            if DC1(1,k)~=0 && UC1(1,k)~=0
                Cell_down=CELL{DC1(1,k)};
                Cell_up=CELL{UC1(1,k)};
                n_u2d=(Cell_down{5}-Cell_up{5})';
                if single(e+n_u2d*V1(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            elseif DC1(1,k)==0 && UC1(1,k)~=0
                Cell_up=CELL{UC1(1,k)};
                n_u2d=(FC{7}-Cell_up{5})';
                if single(e+n_u2d*V1(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            elseif DC1(1,k)~=0 && UC1(1,k)==0
                Cell_down=CELL{DC1(1,k)};
                n_u2d=(Cell_down{5}-FC{7})';
                if single(e+n_u2d*V1(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            else
                error('There is error in FC{12} and FC{13}!');
            end
        end
        % S{1}, S{2} for V2
        DC2=zeros(1,q2);
        UC2=zeros(1,q2);
        for k=1:q2
            if FC{4}*V2(:,k)>=0
                DC2(1,k)=neigh_down(1,1);
                UC2(1,k)=neigh_up(1,1);
            else
                DC2(1,k)=neigh_up(1,1);
                UC2(1,k)=neigh_down(1,1);
            end
        end
        % Check
        % Cell pair
        for k=1:q2
            if length(unique(union([DC2(1,k),UC2(1,k)],[neigh_down(1,1),neigh_up(1,1)])))~=length([neigh_down(1,1),neigh_up(1,1)])
                error('The downwind & upwind cells are not correct!')
            end
        end
        % Direction
        for k=1:q2
            if DC2(1,k)~=0 && UC2(1,k)~=0
                Cell_down=CELL{DC2(1,k)};
                Cell_up=CELL{UC2(1,k)};
                n_u2d=(Cell_down{5}-Cell_up{5})';
                if single(e+n_u2d*V2(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            elseif DC2(1,k)==0 && UC2(1,k)~=0
                Cell_up=CELL{UC2(1,k)};
                n_u2d=(FC{7}-Cell_up{5})';
                if single(e+n_u2d*V2(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            elseif DC2(1,k)~=0 && UC2(1,k)==0
                Cell_down=CELL{DC2(1,k)};
                n_u2d=(Cell_down{5}-FC{7})';
                if single(e+n_u2d*V2(:,k))<single(e)
                    error('The downwind & upwind cells should be switched!');
                end
            else
                error('There is error in FC{12} and FC{13}!');
            end
        end
        % Fill the data
        S_V1{1,1}=DC1;
        S_V1{2,1}=UC1;
        
        S_V2{1,1}=DC2;
        S_V2{2,1}=UC2;
        %% S{4}~S{7}
        %%%%%%% The output of this section C_b, C_j, C_j_down, C_j_up, d_cc, d_nc, dc_down, dc_up can be used multiple times
        if FC{2}~=0 % current face is on boundary
            C_b=FC{7};
            neigh_up=FC{12};
            neigh_down=FC{13};
            cell_in=setxor(0,[neigh_down(1,1),neigh_up(1,1)]);
            if length(cell_in)~=1
                error('There is only one cell attached to the boundary face!');
            end
            P=CELL{cell_in};
            C_j_d=norm_joint(C_b,FC{4}',P{5});
            if in_triangle(C_j_d,P{13},P{14},P{15})
                ;
            else
                error('The generated point is not within the current cell!');
            end
            d_cc=dis(C_b,C_j_d); % distance from the face midpoint
            % to the point that is closest to the centroid of current cell
            d_nc=d_cc; % distance from the face midpoint
            % to the point that is closest to the centroid of neighbor cell
        else  % Internal edge
            C_b=FC{7};
            neigh_up=FC{12};
            neigh_down=FC{13};
            P=CELL{neigh_down(1,1)};
            Q=CELL{neigh_up(1,1)};
            C_j_down=norm_joint(C_b,FC{4}',P{5});
            C_j_up=norm_joint(C_b,FC{4}',Q{5});
            if in_triangle(C_j_down,P{13},P{14},P{15})
                ;
            else
                error('The generated point is not within the current cell!');
            end
            if in_triangle(C_j_up,Q{13},Q{14},Q{15})
                ;
            else
                error('The generated point is not within the current cell!');
            end
            d_c_down=dis(C_b,C_j_down);
            d_c_up=dis(C_b,C_j_up);
        end
        %%%%%%% The output of this section C_b, C_j, C_j_down, C_j_up, d_cc, d_nc, dc_down, dc_up can be used multiple times
        % S{4}~S{7} for V1
        DC1=S_V1{1,1};
        UC1=S_V1{2,1};
        x_down1=zeros(1,q1);
        y_down1=zeros(1,q1);
        x_up1=zeros(1,q1);
        y_up1=zeros(1,q1);
        for k=1:q1
            if DC1(1,k)==neigh_down(1,1) && UC1(1,k)==neigh_up(1,1)
                if FC{2}~=0
                    if DC1(1,k)==0 && UC1(1,k)~=0
                        x_up1(1,k)=C_j_d(1,1);
                        y_up1(1,k)=C_j_d(2,1);
                        x_down1(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_down1(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_down1(1,k)<X2 && x_down1(1,k)>X1) && (y_down1(1,k)<Y2 && y_down1(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    elseif DC1(1,k)~=0 && UC1(1,k)==0
                        x_down1(1,k)=C_j_d(1,1);
                        y_down1(1,k)=C_j_d(2,1);
                        x_up1(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_up1(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_up1(1,k)<X2 && x_up1(1,k)>X1) && (y_up1(1,k)<Y2 && y_up1(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    else
                        error('The face on the boundary contains false neighbor cell info!');
                    end
                else
                    x_down1(1,k)=C_j_down(1,1);
                    y_down1(1,k)=C_j_down(2,1);
                    x_up1(1,k)=C_j_up(1,1);
                    y_up1(1,k)=C_j_up(2,1);
                    if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                    if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                end
            elseif DC1(1,k)==neigh_up(1,1) && UC1(1,k)==neigh_down(1,1)
                if FC{2}~=0
                    if DC1(1,k)==0 && UC1(1,k)~=0
                        x_up1(1,k)=C_j_d(1,1);
                        y_up1(1,k)=C_j_d(2,1);
                        x_down1(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_down1(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_down1(1,k)<X2 && x_down1(1,k)>X1) && (y_down1(1,k)<Y2 && y_down1(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    elseif DC1(1,k)~=0 && UC1(1,k)==0
                        x_down1(1,k)=C_j_d(1,1);
                        y_down1(1,k)=C_j_d(2,1);
                        x_up1(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_up1(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_up1(1,k)<X2 && x_up1(1,k)>X1) && (y_up1(1,k)<Y2 && y_up1(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    else
                        error('The face on the boundary contains false neighbor cell info!');
                    end
                else
                    x_down1(1,k)=C_j_up(1,1);
                    y_down1(1,k)=C_j_up(2,1);
                    x_up1(1,k)=C_j_down(1,1);
                    y_up1(1,k)=C_j_down(2,1);
                    if (x_down1(1,k)<=X1 || x_down1(1,k)>=X2) || (y_down1(1,k)<=Y1 || y_down1(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                    if (x_up1(1,k)<=X1 || x_up1(1,k)>=X2) || (y_up1(1,k)<=Y1 || y_up1(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                end
            else
                error('Two neighbor cells for the current face reside both on one side of the face!');
            end
        end
        % Fill data S{4}~S{7} for V1
        S_V1{4,1}=x_down1;
        S_V1{5,1}=y_down1;
        S_V1{6,1}=x_up1;
        S_V1{7,1}=y_up1;
        
        % S{4}~S{7} for V2
        DC2=S_V2{1,1};
        UC2=S_V2{2,1};
        x_down2=zeros(1,q2);
        y_down2=zeros(1,q2);
        x_up2=zeros(1,q2);
        y_up2=zeros(1,q2);
        for k=1:q2
            if DC2(1,k)==neigh_down(1,1) && UC2(1,k)==neigh_up(1,1)
                if FC{2}~=0
                    if DC2(1,k)==0 && UC2(1,k)~=0
                        x_up2(1,k)=C_j_d(1,1);
                        y_up2(1,k)=C_j_d(2,1);
                        x_down2(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_down2(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_down2(1,k)<X2 && x_down2(1,k)>X1) && (y_down2(1,k)<Y2 && y_down2(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    elseif DC2(1,k)~=0 && UC2(1,k)==0
                        x_down2(1,k)=C_j_d(1,1);
                        y_down2(1,k)=C_j_d(2,1);
                        x_up2(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_up2(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_up2(1,k)<X2 && x_up2(1,k)>X1) && (y_up2(1,k)<Y2 && y_up2(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    else
                        error('The face on the boundary contains false neighbor cell info!');
                    end
                else
                    x_down2(1,k)=C_j_down(1,1);
                    y_down2(1,k)=C_j_down(2,1);
                    x_up2(1,k)=C_j_up(1,1);
                    y_up2(1,k)=C_j_up(2,1);
                    if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                    if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                end
            elseif DC2(1,k)==neigh_up(1,1) && UC2(1,k)==neigh_down(1,1)
                if FC{2}~=0
                    if DC2(1,k)==0 && UC2(1,k)~=0
                        x_up2(1,k)=C_j_d(1,1);
                        y_up2(1,k)=C_j_d(2,1);
                        x_down2(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_down2(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_down2(1,k)<X2 && x_down2(1,k)>X1) && (y_down2(1,k)<Y2 && y_down2(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    elseif DC2(1,k)~=0 && UC2(1,k)==0
                        x_down2(1,k)=C_j_d(1,1);
                        y_down2(1,k)=C_j_d(2,1);
                        x_up2(1,k)=2*C_b(1,1)-C_j_d(1,1);
                        y_up2(1,k)=2*C_b(2,1)-C_j_d(2,1);
                        if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
                            error('The stencil point for one side of the face on boundary should be located inside of the computational domain!');
                        end
                        if (x_up2(1,k)<X2 && x_up2(1,k)>X1) && (y_up2(1,k)<Y2 && y_up2(1,k)>Y1)
                            error('The stencil point for one side of the face on boundary should be located outside of the computational domain!');
                        end
                    else
                        error('The face on the boundary contains false neighbor cell info!');
                    end
                else
                    x_down2(1,k)=C_j_up(1,1);
                    y_down2(1,k)=C_j_up(2,1);
                    x_up2(1,k)=C_j_down(1,1);
                    y_up2(1,k)=C_j_down(2,1);
                    if (x_down2(1,k)<=X1 || x_down2(1,k)>=X2) || (y_down2(1,k)<=Y1 || y_down2(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                    if (x_up2(1,k)<=X1 || x_up2(1,k)>=X2) || (y_up2(1,k)<=Y1 || y_up2(1,k)>=Y2)
                        error('The stencil point for any side of the interior face should be located inside of the computational domain!');
                    end
                end
            else
                error('Two neighbor cells for the current face reside both on one side of the face!');
            end
        end
        % Fill data S{4}~S{7} for V2
        S_V2{4,1}=x_down2;
        S_V2{5,1}=y_down2;
        S_V2{6,1}=x_up2;
        S_V2{7,1}=y_up2;
        
        %% S{13}, S{15}
        % S{13}, S{15} for V1
        x_down1=S_V1{4,1};
        y_down1=S_V1{5,1};
        x_up1=S_V1{6,1};
        y_up1=S_V1{7,1};
        dis_u2d=zeros(1,q1);
        One_over_dis_u2d=zeros(1,q1);
        for k=1:q1
            dis_u2d(1,k)=dis([x_down1(1,k);y_down1(1,k)],[x_up1(1,k);y_up1(1,k)]);
            if FC{2}~=0
                if single(e+dis_u2d(1,k))~=single(e+2*d_nc)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            else
                if single(e+dis_u2d(1,k))~=single(e+d_c_down+d_c_up)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            end
            One_over_dis_u2d(1,k)=1/dis_u2d(1,k);
        end
        % Fill the data S{13}, S{15} for V1
        S_V1{13,1}=dis_u2d;
        S_V1{15,1}=One_over_dis_u2d;
        
        % S{13}, S{15} for V2
        x_down2=S_V2{4,1};
        y_down2=S_V2{5,1};
        x_up2=S_V2{6,1};
        y_up2=S_V2{7,1};
        dis_u2d=zeros(1,q2);
        One_over_dis_u2d=zeros(1,q2);
        for k=1:q2
            dis_u2d(1,k)=dis([x_down2(1,k);y_down2(1,k)],[x_up2(1,k);y_up2(1,k)]);
            if FC{2}~=0
                if single(e+dis_u2d(1,k))~=single(e+2*d_nc)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            else
                if single(e+dis_u2d(1,k))~=single(e+d_c_down+d_c_up)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            end
            One_over_dis_u2d(1,k)=1/dis_u2d(1,k);
        end
        % Fill the data S{13}, S{15} for V2
        S_V2{13,1}=dis_u2d;
        S_V2{15,1}=One_over_dis_u2d;
        
        %% S{17}, S{19}
        % S{17}, S{19} for V1
        x_down1=S_V1{4,1};
        y_down1=S_V1{5,1};
        x_up1=S_V1{6,1};
        y_up1=S_V1{7,1};
        u2d_over_u=zeros(1,q1);
        One_u2d_over_u=zeros(1,q1);
        for k=1:q1
            u2d_over_u(1,k)=dis([x_down1(1,k);y_down1(1,k)],[x_up1(1,k);y_up1(1,k)])/dis(C_b,[x_up1(1,k);y_up1(1,k)]);
            if FC{2}~=0
                if single(e+u2d_over_u(1,k))~=single(e+2)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            else
                if single(e+u2d_over_u(1,k))~=single(e+(d_c_down+d_c_up)/dis(C_b,[x_up1(1,k);y_up1(1,k)]))
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            end
            One_u2d_over_u(1,k)=1/u2d_over_u(1,k);
        end
        % Fill the data S{13}, S{15} for V1
        S_V1{17,1}=u2d_over_u;
        S_V1{19,1}=One_u2d_over_u;
        
        % S{17}, S{19} for V2
        x_down2=S_V2{4,1};
        y_down2=S_V2{5,1};
        x_up2=S_V2{6,1};
        y_up2=S_V2{7,1};
        u2d_over_u=zeros(1,q2);
        One_u2d_over_u=zeros(1,q2);
        for k=1:q2
            u2d_over_u(1,k)=dis([x_down2(1,k);y_down2(1,k)],[x_up2(1,k);y_up2(1,k)])/dis(C_b,[x_up2(1,k);y_up2(1,k)]);
            if FC{2}~=0
                if single(e+u2d_over_u(1,k))~=single(e+2)
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            else
                if single(e+u2d_over_u(1,k))~=single(e+(d_c_down+d_c_up)/dis(C_b,[x_up2(1,k);y_up2(1,k)]))
                    error('The distance from upwind to downwind cell on the for the boundary face is incorrect!');
                end
            end
            One_u2d_over_u(1,k)=1/u2d_over_u(1,k);
        end
        % Fill the data S{13}, S{15} for V1
        S_V2{17,1}=u2d_over_u;
        S_V2{19,1}=One_u2d_over_u;
        %% FC{18}
        fc18=zeros(2,3);
        if FC{2}==0
            fc18(:,1)=C_j_down;
            fc18(:,2)=C_j_up;
        else
            if neigh_up(1,1)==0 && neigh_down(1,1)~=0
                fc18(:,1)=C_j_d;
            elseif neigh_up(1,1)~=0 && neigh_down(1,1)==0
                fc18(:,2)=C_j_d;
            else
                error('There is error in FC{12} and FC{13}!');
            end
        end
        FC{18}=fc18;
        %% S{10}, S{11}
        neigh_up=FC{12};
        neigh_down=FC{13};
        ZONE_UP=-1;
        ZONE_DOWN=-1;
        if FC{2}~=0
            if neigh_down(1,1)==0 && neigh_up(1,1)~=0
                Cell_up=CELL{neigh_up(1,1)};
                nd1_up=Cell_up{13};
                nd2_up=Cell_up{14};
                nd3_up=Cell_up{15};
                in_zone_one=in_triangle(C_j_d,Cell_up{5},nd1_up,nd2_up);
                in_zone_two=in_triangle(C_j_d,Cell_up{5},nd2_up,nd3_up);
                in_zone_three=in_triangle(C_j_d,Cell_up{5},nd3_up,nd1_up);
                if in_zone_one
                    ZONE_UP=1; % Located in zone 1
                elseif in_zone_two
                    ZONE_UP=2; % Located in zone 2
                elseif in_zone_three
                    ZONE_UP=3; % Located in zone 3
                else
                    on_edge_one=on_edge(C_j_d,Cell_up{5},nd1_up);
                    on_edge_two=on_edge(C_j_d,Cell_up{5},nd2_up);
                    on_edge_three=on_edge(C_j_d,Cell_up{5},nd3_up);
                    if on_edge_one && on_edge_two && on_edge_three
                        ZONE_UP=4;
                    elseif on_edge_one
                        ZONE_UP=1;
                    elseif on_edge_two
                        ZONE_UP=2;
                    elseif on_edge_three
                        ZONE_UP=3;
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
                ZONE_DOWN=4;
            elseif neigh_down(1,1)~=0 && neigh_up(1,1)==0
                Cell_down=CELL{neigh_down(1,1)};
                nd1_down=Cell_down{13};
                nd2_down=Cell_down{14};
                nd3_down=Cell_down{15};
                in_zone_one=in_triangle(C_j_d,Cell_down{5},nd1_down,nd2_down);
                in_zone_two=in_triangle(C_j_d,Cell_down{5},nd2_down,nd3_down);
                in_zone_three=in_triangle(C_j_d,Cell_down{5},nd3_down,nd1_down);
                if in_zone_one
                    ZONE_DOWN=1; % Located in zone 1
                elseif in_zone_two
                    ZONE_DOWN=2; % Located in zone 2
                elseif in_zone_three
                    ZONE_DOWN=3; % Located in zone 3
                else
                    on_edge_one=on_edge(C_j_d,Cell_down{5},nd1_down);
                    on_edge_two=on_edge(C_j_d,Cell_down{5},nd2_down);
                    on_edge_three=on_edge(C_j_d,Cell_down{5},nd3_down);
                    if on_edge_one && on_edge_two && on_edge_three
                        ZONE_DOWN=4;
                    elseif on_edge_one
                        ZONE_DOWN=1;
                    elseif on_edge_two
                        ZONE_DOWN=2;
                    elseif on_edge_three
                        ZONE_DOWN=3;
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
                ZONE_UP=4;
            else
                error('The face on the boundary contains false neighbor cell info!');
            end
        else
            if FC{10}~=0 || FC{11}~=0
                Cell_down=CELL{neigh_down(1,1)};
                Cell_up=CELL{neigh_up(1,1)};
                nd1_down=Cell_down{13};
                nd2_down=Cell_down{14};
                nd3_down=Cell_down{15};
                nd1_up=Cell_up{13};
                nd2_up=Cell_up{14};
                nd3_up=Cell_up{15};
                % down cell
                in_zone_one=in_triangle(C_j_down,Cell_down{5},nd1_down,nd2_down);
                in_zone_two=in_triangle(C_j_down,Cell_down{5},nd2_down,nd3_down);
                in_zone_three=in_triangle(C_j_down,Cell_down{5},nd3_down,nd1_down);
                if in_zone_one
                    ZONE_DOWN=1; % Located in zone 1
                elseif in_zone_two
                    ZONE_DOWN=2; % Located in zone 2
                elseif in_zone_three
                    ZONE_DOWN=3; % Located in zone 3
                else
                    on_edge_one=on_edge(C_j_down,Cell_down{5},nd1_down);
                    on_edge_two=on_edge(C_j_down,Cell_down{5},nd2_down);
                    on_edge_three=on_edge(C_j_down,Cell_down{5},nd3_down);
                    if on_edge_one && on_edge_two && on_edge_three
                        ZONE_DOWN=4;
                    elseif on_edge_one
                        ZONE_DOWN=1;
                    elseif on_edge_two
                        ZONE_DOWN=2;
                    elseif on_edge_three
                        ZONE_DOWN=3;
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
                % up cell
                in_zone_one=in_triangle(C_j_up,Cell_up{5},nd1_up,nd2_up);
                in_zone_two=in_triangle(C_j_up,Cell_up{5},nd2_up,nd3_up);
                in_zone_three=in_triangle(C_j_up,Cell_up{5},nd3_up,nd1_up);
                if in_zone_one
                    ZONE_UP=1; % Located in zone 1
                elseif in_zone_two
                    ZONE_UP=2; % Located in zone 2
                elseif in_zone_three
                    ZONE_UP=3; % Located in zone 3
                else
                    on_edge_one=on_edge(C_j_up,Cell_up{5},nd1_up);
                    on_edge_two=on_edge(C_j_up,Cell_up{5},nd2_up);
                    on_edge_three=on_edge(C_j_up,Cell_up{5},nd3_up);
                    if on_edge_one && on_edge_two && on_edge_three
                        ZONE_UP=4;
                    elseif on_edge_one
                        ZONE_UP=1;
                    elseif on_edge_two
                        ZONE_UP=2;
                    elseif on_edge_three
                        ZONE_UP=3;
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
            end
        end
        %     %  S{10}, S{11} for V1
        %     DC1=S_V1{1,1};
        %     UC1=S_V1{2,1};
        %     zone_up1=zeros(1,q1);
        %     zone_down1=zeros(1,q1);
        %     for k=1:q1
        %         if DC1(1,k)==neigh_down(1,1) && UC1(1,k)==neigh_up(1,1)
        %             zone_down1(1,k)=ZONE_DOWN;
        %             zone_up1(1,k)=ZONE_UP;
        %         elseif DC1(1,k)==neigh_up(1,1) && UC1(1,k)==neigh_down(1,1)
        %             zone_down1(1,k)=ZONE_UP;
        %             zone_up1(1,k)=ZONE_DOWN;
        %         else
        %             error('Two neighbor cells for the current face reside both on one side of the face!');
        %         end
        %     end
        %     %  Fill the data S{10}, S{11} for V1
        %     S_V1{10,1}=zone_down1;
        %     S_V1{11,1}=zone_up1;
        %
        %     %  S{10}, S{11} for V2
        %     DC2=S_V2{1,1};
        %     UC2=S_V2{2,1};
        %     zone_up2=zeros(1,q2);
        %     zone_down2=zeros(1,q2);
        %     for k=1:q2
        %         if DC2(1,k)==neigh_down(1,1) && UC2(1,k)==neigh_up(1,1)
        %             zone_down2(1,k)=ZONE_DOWN;
        %             zone_up2(1,k)=ZONE_UP;
        %         elseif DC2(1,k)==neigh_up(1,1) && UC2(1,k)==neigh_down(1,1)
        %             zone_down2(1,k)=ZONE_UP;
        %             zone_up2(1,k)=ZONE_DOWN;
        %         else
        %             error('Two neighbor cells for the current face reside both on one side of the face!');
        %         end
        %     end
        %     %  Fill the data S{10}, S{11} for V2
        %     S_V2{10,1}=zone_down2;
        %     S_V2{11,1}=zone_up2;
        %% Fill in FC{16}, FC{17}
        FC{16}=S_V1;
        FC{17}=S_V2;
        %% FC{19}
        fc19=zeros(1,3);
        fc19(:,1)=ZONE_DOWN;
        fc19(:,2)=ZONE_UP;
        FC{19}=fc19;
        %% FC{20}
        fc20=zeros(4,3);
        if FC{10}==0 && FC{11}==0 % Interior face, and none of the end nodes is on boundary
            UP=FC{12};
            DOWN=FC{13};
            UP_Cell_number=UP(1,1);
            DOWN_Cell_number=DOWN(1,1);
            Cell_up_fixed=CELL{UP_Cell_number};
            C_up_fixed=Cell_up_fixed{5};
            Cell_down_fixed=CELL{DOWN_Cell_number};
            C_down_fixed=Cell_down_fixed{5};
            ND1=NODE{FC{8}};
            ND2=NODE{FC{9}};
            Cell_union=union(ND1{5},ND2{5});
            Cell_union_UP=setxor(UP_Cell_number,Cell_union);
            Cell_union_DOWN=setxor(DOWN_Cell_number,Cell_union);
            L=length(Cell_union_UP);
            L_cell_union_down=length(Cell_union_DOWN);
            Pair_cell_union_up=zeros(2,L*(L-1)/2);
            Pair_cell_union=zeros(2,L_cell_union_down*(L_cell_union_down-1)/2);
            S=FC{18};
            C_down=S(:,1);
            C_up=S(:,2);
            %% Find the triangle that circles the stencil point in the downstream cell
            if single(e+norm(C_down_fixed-C_down))==single(e) % The stencil point is at the centroid
                fc20(1,1)=0;
                fc20(2,1)=DOWN_Cell_number;
                fc20(3,1)=DOWN_Cell_number;
                fc20(4,1)=DOWN_Cell_number;
            else
                %% Find all possible paired centroids that could be used to form triangles
                a=0;
                for k=1:length(Cell_union_DOWN)
                    if length(Cell_union_DOWN)==2
                        a=a+1;
                        Pair_cell_union(:,a)=[Cell_union_DOWN(1);Cell_union_DOWN(2)];
                        break;
                    else
                        for n=2:length(Cell_union_DOWN)
                            a=a+1;
                            Pair_cell_union(:,a)=[Cell_union_DOWN(1);Cell_union_DOWN(n)];
                        end
                        Cell_union_DOWN=setxor(Cell_union_DOWN(1),Cell_union_DOWN);
                    end
                end
                %% Narrow down the previous possibility ruling out the triangles that don't
                % circle the stencil point in the downstream cell
                a=0;
                in_tri=0;
                for k=1:a
                    Cell_pair_1=CELL{Pair_cell_union(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union(2,k)};
                    C_down_1=Cell_pair_1{5};
                    C_down_2=Cell_pair_2{5};
                    if in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                        a=a+1;
                        in_tri(a)=k;
                    end
                end
                if a==0
                    error('No triangle is found that contains the stencil point!');
                end
                
                Pair_cell_union_down_new=zeros(2,a);
                for n=1:a
                    Pair_cell_union_down_new(:,n)=Pair_cell_union(:,in_tri(n));
                end
                % Check
                for k=1:a
                    Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
                    C_down_1=Cell_pair_1{5};
                    C_down_2=Cell_pair_2{5};
                    if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                        error(':');
                    end
                end
                %% Find the pair that has the shortest distance to the stencil point
                Dis=zeros(1,a);
                for k=1:a
                    Cell_pair_1=CELL{Pair_cell_union_down_new(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union_down_new(2,k)};
                    C_down_1=Cell_pair_1{5};
                    C_down_2=Cell_pair_2{5};
                    Dis(1,k)=(dis(C_down,C_down_1)+dis(C_down,C_down_2))/2;
                end
                Dis_min=min(Dis);
                for k=1:a
                    if single(e+Dis_min)==single(e+Dis(1,k))
                        break;
                    end
                end
                Cell_pair_down_found=Pair_cell_union_down_new(:,k);
                % Check
                Cell_pair_1=CELL{Cell_pair_down_found(1,1)};
                Cell_pair_2=CELL{Cell_pair_down_found(2,1)};
                C_down_1=Cell_pair_1{5};
                C_down_2=Cell_pair_2{5};
                if ~in_triangle(C_down,C_down_fixed,C_down_1,C_down_2)
                    error(':');
                end
                fc20(1,1)=0;
                fc20(2,1)=DOWN_Cell_number;
                fc20(3,1)=Cell_pair_down_found(1,1);
                fc20(4,1)=Cell_pair_down_found(2,1);
            end
            %% Find the triangle that circles the stencil point in the upstream cell
            if single(e+norm(C_up_fixed-C_up))==single(e) % The stencil point is at the cnetroid
                fc20(1,2)=0;
                fc20(2,2)=UP_Cell_number;
                fc20(3,2)=UP_Cell_number;
                fc20(4,2)=UP_Cell_number;
            else
                %% Find all possible paired centroids that could be used to form triangles
                a_up=0;
                for k=1:length(Cell_union_UP)
                    if length(Cell_union_UP)==2
                        a_up=a_up+1;
                        Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(2)];
                        break;
                    else
                        for n=2:length(Cell_union_UP)
                            a_up=a_up+1;
                            Pair_cell_union_up(:,a_up)=[Cell_union_UP(1);Cell_union_UP(n)];
                        end
                        Cell_union_UP=setxor(Cell_union_UP(1),Cell_union_UP);
                    end
                end
                %% Narrow down the previous possibility ruling out the triangles that don't
                % circle the stencil point in the upstream cell
                a=0;
                in_tri_pair=0;
                for k=1:a_up
                    Cell_pair_1=CELL{Pair_cell_union_up(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union_up(2,k)};
                    C_up_1=Cell_pair_1{5};
                    C_up_2=Cell_pair_2{5};
                    if in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                        a=a+1;
                        in_tri_pair(a)=k;
                    end
                end
                if a==0
                    error('No triangle is found that contains the stencil point!');
                end
                
                Pair_cell_union_up_new=zeros(2,a);
                for n=1:a
                    Pair_cell_union_up_new(:,n)=Pair_cell_union_up(:,in_tri_pair(n));
                end
                % Check
                for k=1:a
                    Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
                    C_up_1=Cell_pair_1{5};
                    C_up_2=Cell_pair_2{5};
                    if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                        error(':');
                    end
                end
                %% Find the pair that has the shortest distance to the stencil point
                Dis_up=zeros(1,a);
                for k=1:a
                    Cell_pair_1=CELL{Pair_cell_union_up_new(1,k)};
                    Cell_pair_2=CELL{Pair_cell_union_up_new(2,k)};
                    C_up_1=Cell_pair_1{5};
                    C_up_2=Cell_pair_2{5};
                    Dis_up(1,k)=(dis(C_up,C_up_1)+dis(C_up,C_up_2))/2;
                end
                Dis_up_min=min(Dis_up);
                for k=1:a
                    if single(e+Dis_up_min)==single(e+Dis_up(1,k))
                        break;
                    end
                end
                Cell_pair_up_found=Pair_cell_union_up_new(:,k);
                % Check
                Cell_pair_1=CELL{Cell_pair_up_found(1,1)};
                Cell_pair_2=CELL{Cell_pair_up_found(2,1)};
                C_up_1=Cell_pair_1{5};
                C_up_2=Cell_pair_2{5};
                if ~in_triangle(C_up,C_up_fixed,C_up_1,C_up_2)
                    error(':');
                end
                fc20(1,2)=0;
                fc20(2,2)=UP_Cell_number;
                fc20(3,2)=Cell_pair_up_found(1,1);
                fc20(4,2)=Cell_pair_up_found(2,1);
            end
        elseif (FC{10}~=0 && FC{11}==0) || (FC{10}==0 && FC{11}~=0) % Interior face, one of the end nodes is on boundary, the other is not
            % Check
            if FC{2}~=0
                error('The current face should be interior!');
            end
            UP=FC{12};
            DOWN=FC{13};
            UP_Cell_number=UP(1,1);
            DOWN_Cell_number=DOWN(1,1);
            Cell_up_fixed=CELL{UP_Cell_number};
            C_up_fixed=Cell_up_fixed{5};
            Cell_down_fixed=CELL{DOWN_Cell_number};
            C_down_fixed=Cell_down_fixed{5};
            ND1=NODE{FC{8}};
            ND2=NODE{FC{9}};
            if FC{10}~=0 && FC{11}==0
                C_nd=ND1{3};
                nd=ND1{1};
                Cell_union=union(ND1{5},ND2{5}); % Chould be the node on boundary
            elseif FC{10}==0 && FC{11}~=0
                C_nd=ND2{3};
                nd=ND2{1};
                Cell_union=union(ND1{5},ND2{5}); % Chould be the node on boundary
            else
                error('There must be one node on the boundary and the other is not!');
            end
            Third_cell_up=setxor(UP_Cell_number,Cell_union);
            Third_cell_down=setxor(DOWN_Cell_number,Cell_union);
            S=FC{18};
            C_down=S(:,1);
            C_up=S(:,2);
            %% Find the triangle that circles the stencil point in the downstream cell
            if single(e+norm(C_down_fixed-C_down))==single(e) % The stencil point is at the cnetroid
                fc20(1,1)=0;
                fc20(2,1)=DOWN_Cell_number;
                fc20(3,1)=DOWN_Cell_number;
                fc20(4,1)=DOWN_Cell_number;
            else
                %% Find all possible paired centroids that could be used to form triangles
                
                %% Narrow down the previous possibility ruling out the triangles that don't
                % circle the stencil point in the downstream cell
                a=length(Third_cell_down);
                a=0;
                in_tri=0;
                for k=1:a
                    Cell_third_pool=CELL{Third_cell_down(k)};
                    C_third_cell=Cell_third_pool{5};
                    if in_triangle(C_down,C_down_fixed,C_nd,C_third_cell)
                        a=a+1;
                        in_tri(a)=k;
                    end
                end
                if a==0
                    if Cell_down_fixed{10}==0 && Cell_down_fixed{11}~=0 && Cell_down_fixed{12}~=0
                        nd1=Cell_down_fixed{8};
                        nd2=Cell_down_fixed{9};
                        C_nd1=Cell_down_fixed{14};
                        C_nd2=Cell_down_fixed{15};
                    elseif Cell_down_fixed{10}~=0 && Cell_down_fixed{11}==0 && Cell_down_fixed{12}~=0
                        nd1=Cell_down_fixed{7};
                        nd2=Cell_down_fixed{9};
                        C_nd1=Cell_down_fixed{13};
                        C_nd2=Cell_down_fixed{15};
                    elseif Cell_down_fixed{10}~=0 && Cell_down_fixed{11}~=0 && Cell_down_fixed{12}==0
                        nd1=Cell_down_fixed{7};
                        nd2=Cell_down_fixed{8};
                        C_nd1=Cell_down_fixed{13};
                        C_nd2=Cell_down_fixed{14};
                    else
                        error('The down cell should be attached to boundary!');
                    end
                    if ~in_triangle(C_down,C_down_fixed,C_nd1,C_nd2)
                        error('The stencil point should be located in the one of the ZONE of the cell attached to boundary!');
                    end
                    fc20(1,1)=2;
                    fc20(2,1)=nd1;
                    fc20(3,1)=nd2;
                    fc20(4,1)=DOWN_Cell_number;
                else
                    Third_cell_down_new=zeros(1,a);
                    for n=1:a
                        Third_cell_down_new(1,n)=Third_cell_down(in_tri(n));
                    end
                    % Check
                    for k=1:a
                        Cell_third_pool=CELL{Third_cell_down_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        if ~in_triangle(C_down,C_down_fixed,C_nd,C_third_cell)
                            error(':');
                        end
                    end
                    %% Find the pair that has the shortest distance to the stencil point
                    Dis=zeros(1,a);
                    for k=1:a
                        Cell_third_pool=CELL{Third_cell_down_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        Dis(1,k)=dis(C_down,C_third_cell);
                    end
                    Dis_min=min(Dis);
                    for k=1:a
                        if single(e+Dis_min)==single(e+Dis(1,k))
                            break;
                        end
                    end
                    Third_cell_found=Third_cell_down_new(1,k);
                    % Check
                    Cell_third_pool=CELL{Third_cell_found};
                    C_third_cell=Cell_third_pool{5};
                    if ~in_triangle(C_down,C_down_fixed,C_nd,C_third_cell)
                        error(':');
                    end
                    fc20(1,1)=1;
                    fc20(2,1)=nd;
                    fc20(3,1)=DOWN_Cell_number;
                    fc20(4,1)=Third_cell_found;
                end
            end
            %% Find the triangle that circles the stencil point in the upstream cell
            if single(e+norm(C_up_fixed-C_up))==single(e) % The stencil point is at the cnetroid
                fc20(1,2)=0;
                fc20(2,2)=UP_Cell_number;
                fc20(3,2)=UP_Cell_number;
                fc20(4,2)=UP_Cell_number;
            else
                %% Find all possible paired centroids that could be used to form triangles
                
                %% Narrow down the previous possibility ruling out the triangles that don't
                % circle the stencil point in the downstream cell
                a_up=length(Third_cell_up);
                a=0;
                in_tri_pair=0;
                for k=1:a_up
                    Cell_third_pool=CELL{Third_cell_up(k)};
                    C_third_cell_up=Cell_third_pool{5};
                    if in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
                        a=a+1;
                        in_tri_pair(a)=k;
                    end
                end
                if a==0
                    if Cell_up_fixed{10}==0 && Cell_up_fixed{11}~=0 && Cell_up_fixed{12}~=0
                        nd1=Cell_up_fixed{8};
                        nd2=Cell_up_fixed{9};
                        C_nd1=Cell_up_fixed{14};
                        C_nd2=Cell_up_fixed{15};
                    elseif Cell_up_fixed{10}~=0 && Cell_up_fixed{11}==0 && Cell_up_fixed{12}~=0
                        nd1=Cell_up_fixed{7};
                        nd2=Cell_up_fixed{9};
                        C_nd1=Cell_up_fixed{13};
                        C_nd2=Cell_up_fixed{15};
                    elseif Cell_up_fixed{10}~=0 && Cell_up_fixed{11}~=0 && Cell_up_fixed{12}==0
                        nd1=Cell_up_fixed{7};
                        nd2=Cell_up_fixed{8};
                        C_nd1=Cell_up_fixed{13};
                        C_nd2=Cell_up_fixed{14};
                    else
                        error('The up cell should be attached to boundary!');
                    end
                    if ~in_triangle(C_up,C_up_fixed,C_nd1,C_nd2)
                        error('The stencil point should be located in the one of the ZONE of the cell attached to boundary!');
                    end
                    fc20(1,2)=2;
                    fc20(2,2)=nd1;
                    fc20(3,2)=nd2;
                    fc20(4,2)=UP_Cell_number;
                else
                    Third_cell_up_new=zeros(1,a);
                    for n=1:a
                        Third_cell_up_new(1,n)=Third_cell_up(in_tri_pair(n));
                    end
                    % Check
                    for k=1:a
                        Cell_third_pool=CELL{Third_cell_up_new(1,k)};
                        C_third_cell_up=Cell_third_pool{5};
                        if ~in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
                            error(':');
                        end
                    end
                    %% Find the pair that has the shortest distance to the stencil point
                    Dis_up=zeros(1,a);
                    for k=1:a
                        Cell_third_pool=CELL{Third_cell_up_new(1,k)};
                        C_third_cell_up=Cell_third_pool{5};
                        Dis_up(1,k)=dis(C_down,C_third_cell_up);
                    end
                    Dis_up_min=min(Dis_up);
                    for k=1:a
                        if single(e+Dis_up_min)==single(e+Dis_up(1,k))
                            break;
                        end
                    end
                    Third_cell_up_found=Third_cell_up_new(1,k);
                    % Check
                    Cell_third_pool=CELL{Third_cell_up_found};
                    C_third_cell_up=Cell_third_pool{5};
                    if ~in_triangle(C_up,C_up_fixed,C_nd,C_third_cell_up)
                        error(':');
                    end
                    fc20(1,2)=1;
                    fc20(2,2)=nd;
                    fc20(3,2)=UP_Cell_number;
                    fc20(4,2)=Third_cell_up_found;
                end
            end
        elseif FC{10}~=0 && FC{11}~=0 % Boundary face
            % Check
            if FC{2}==0
                error('The current face should be interior!');
            end
            fc20_bc=zeros(4,1);
            UP=FC{12};
            DOWN=FC{13};
            UP_Cell_number=UP(1,1);
            DOWN_Cell_number=DOWN(1,1);
            S=FC{18};
            if UP_Cell_number==0 && DOWN_Cell_number~=0
                Cell_number=DOWN_Cell_number;
                C_S=S(:,1);
            elseif UP_Cell_number~=0 && DOWN_Cell_number==0
                Cell_number=UP_Cell_number;
                C_S=S(:,2);
            else
                error('boundary face should have only one neighbor!');
            end
            Cell_fixed=CELL{Cell_number};
            C_fixed=Cell_fixed{5};
            ND1=NODE{FC{8}};
            ND2=NODE{FC{9}};
            ND_3rd=NODE{setxor(union(FC{8},FC{9}),union(Cell_fixed{7},union(Cell_fixed{8},Cell_fixed{9})))};
            
            C_nd1=ND1{3};
            nd1=ND1{1};
            
            C_nd2=ND2{3};
            nd2=ND2{1};
            
            Cell_union=ND_3rd{5};
            if single(e+norm(C_fixed-C_S))==single(e)
                fc20_bc(1,1)=0;
                fc20_bc(2,1)=Cell_number;
                fc20_bc(3,1)=Cell_number;
                fc20_bc(4,1)=Cell_number;
            elseif in_triangle(C_S,C_fixed,C_nd1,C_nd2)
                fc20_bc(1,1)=2;
                fc20_bc(2,1)=FC{8};
                fc20_bc(3,1)=FC{9};
                fc20_bc(4,1)=Cell_number;
            else
                Third_cell=setxor(Cell_number,Cell_union);
                b=length(Third_cell);
                % Using the first end node
                a1=0;
                in_tri_pair_1=0;
                for k=1:b
                    Cell_third_pool=CELL{Third_cell(k)};
                    C_third_cell=Cell_third_pool{5};
                    if in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
                        a1=a1+1;
                        in_tri_pair_1(a1)=k;
                    end
                end
                % Using the second end node
                a2=0;
                in_tri_pair_2=0;
                for k=1:b
                    Cell_third_pool=CELL{Third_cell(k)};
                    C_third_cell=Cell_third_pool{5};
                    if in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
                        a2=a2+1;
                        in_tri_pair_2(a2)=k;
                    end
                end
                
                if a1==0 && a2==0
                    error('There must be a triangle that can encircle the stendil point!');
                elseif a1==0 && a2~=0
                    Third_cell_new=zeros(1,a2);
                    for n=1:a2
                        Third_cell_new(1,n)=Third_cell(in_tri_pair_2(n));
                    end
                    % Check
                    for k=1:a2
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
                            error(':');
                        end
                    end
                    %% Find the pair that has the shortest distance to the stencil point
                    Dis=zeros(1,a2);
                    for k=1:a2
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        Dis(1,k)=dis(C_S,C_third_cell);
                    end
                    Dis_min=min(Dis);
                    for k=1:a2
                        if single(e+Dis_min)==single(e+Dis(1,k))
                            break;
                        end
                    end
                    Third_cell_found=Third_cell_new(1,k);
                    % Check
                    Cell_third_pool=CELL{Third_cell_found};
                    C_third_cell=Cell_third_pool{5};
                    if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
                        error(':');
                    end
                    fc20_bc(1,1)=1;
                    fc20_bc(2,1)=nd2;
                    fc20_bc(3,1)=Cell_number;
                    fc20_bc(4,1)=Third_cell_found;
                elseif a1~=0 && a2==0
                    Third_cell_new=zeros(1,a1);
                    for n=1:a1
                        Third_cell_new(1,n)=Third_cell(in_tri_pair_1(n));
                    end
                    % Check
                    for k=1:a1
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
                            error(':');
                        end
                    end
                    %% Find the pair that has the shortest distance to the stencil point
                    Dis=zeros(1,a1);
                    for k=1:a1
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        Dis(1,k)=dis(C_S,C_third_cell);
                    end
                    Dis_min=min(Dis);
                    for k=1:a1
                        if single(e+Dis_min)==single(e+Dis(1,k))
                            break;
                        end
                    end
                    Third_cell_found=Third_cell_new(1,k);
                    % Check
                    Cell_third_pool=CELL{Third_cell_found};
                    C_third_cell=Cell_third_pool{5};
                    if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
                        error(':');
                    end
                    fc20_bc(1,1)=1;
                    fc20_bc(2,1)=nd1;
                    fc20_bc(3,1)=Cell_number;
                    fc20_bc(4,1)=Third_cell_found;
                elseif a1~=0 && a2~=0
                    Third_cell_new=zeros(1,a1+a2);
                    for n=1:a1
                        Third_cell_new(1,n)=Third_cell(in_tri_pair_1(n));
                    end
                    for n=1:a2
                        Third_cell_new(1,a1+n)=Third_cell(in_tri_pair_2(n));
                    end
                    % Check
                    for k=1:a1+a2
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        if k<=a1
                            if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
                                error(':');
                            end
                        else
                            if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
                                error(':');
                            end
                        end
                    end
                    %% Find the pair that has the shortest distance to the stencil point
                    Dis=zeros(1,a1+a2);
                    for k=1:a1+a2
                        Cell_third_pool=CELL{Third_cell_new(1,k)};
                        C_third_cell=Cell_third_pool{5};
                        Dis(1,k)=dis(C_S,C_third_cell);
                    end
                    Dis_min=min(Dis);
                    for k=1:a1+a2
                        if single(e+Dis_min)==single(e+Dis(1,k))
                            break;
                        end
                    end
                    Third_cell_found=Third_cell_new(1,k);
                    % Check
                    Cell_third_pool=CELL{Third_cell_found};
                    C_third_cell=Cell_third_pool{5};
                    if k<=a1
                        if ~in_triangle(C_S,C_fixed,C_nd1,C_third_cell)
                            error(':');
                        end
                        nd=nd1;
                    else
                        if ~in_triangle(C_S,C_fixed,C_nd2,C_third_cell)
                            error(':');
                        end
                        nd=nd2;
                    end
                    fc20_bc(1,1)=1;
                    fc20_bc(2,1)=nd;
                    fc20_bc(3,1)=Cell_number;
                    fc20_bc(4,1)=Third_cell_found;
                else
                    ;
                end
            end
            if UP_Cell_number==0 && DOWN_Cell_number~=0
                fc20(:,1)=fc20_bc;
            elseif UP_Cell_number~=0 && DOWN_Cell_number==0
                fc20(:,2)=fc20_bc;
            else
                error('boundary face should have only one neighbor!');
            end
            
        else
            error('Boundary identifier of end nodes on the current face is incorrect!');
        end
        FC{20}=fc20;
        %% FC{21}, FC{22}
        fc21=zeros(q1,3);
        fc22=zeros(q2,3);
        if FC{2}==0 % Interior face
            UP=FC{12};
            DOWN=FC{13};
            UP_Cell_number=UP(1,1);
            DOWN_Cell_number=DOWN(1,1);
            P_up=CELL{UP_Cell_number};
            P_down=CELL{DOWN_Cell_number};
            S=FC{18};
            Coord_down=S(:,1);
            Coord_up=S(:,2);
            v_up=Coord_up-P_up{5};
            v_down=Coord_down-P_down{5};
            %% Downwind stencil point
            %% V1
            for s=1:q1
                if V1(:,s)'*v_down>=0
                    fc21(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
            %% V2
            for s=1:q2
                if V2(:,s)'*v_down>=0
                    fc22(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
            %% upwind stencil point
            %% V1
            for s=1:q1
                if V1(:,s)'*v_up>=0
                    fc21(s,2)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
            %% V2
            for s=1:q2
                if V2(:,s)'*v_up>=0
                    fc22(s,2)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
        else % Boundary face
            fc21_bc=zeros(q1,1);
            fc22_bc=zeros(q2,1);
            UP=FC{12};
            DOWN=FC{13};
            UP_Cell_number=UP(1,1);
            DOWN_Cell_number=DOWN(1,1);
            S=FC{18};
            if UP_Cell_number==0 && DOWN_Cell_number~=0
                P=CELL{DOWN_Cell_number};
                Coord=S(:,1);
            elseif UP_Cell_number~=0 && DOWN_Cell_number==0
                P=CELL{UP_Cell_number};
                Coord=S(:,2);
            else
                error('boundary face should have only one neighbor!');
            end
            v=Coord-P{5};
            %% V1
            for s=1:q1
                if V1(:,s)'*v>=0
                    fc21_bc(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
            %% V2
            for s=1:q2
                if V2(:,s)'*v>=0
                    fc22_bc(s,1)=1; % 1 means using the centroid value; 0 means using the 2nd-order mapping value
                end
            end
            if UP_Cell_number==0 && DOWN_Cell_number~=0
                fc21(:,1)=fc21_bc;
                fc22(:,1)=fc22_bc;
            elseif UP_Cell_number~=0 && DOWN_Cell_number==0
                fc21(:,2)=fc21_bc;
                fc22(:,2)=fc22_bc;
            else
                error('boundary face should have only one neighbor!');
            end
        end
        FC{21}=fc21;
        FC{22}=fc22;
        %% Final data filling
        FACE{l}=FC;
    end