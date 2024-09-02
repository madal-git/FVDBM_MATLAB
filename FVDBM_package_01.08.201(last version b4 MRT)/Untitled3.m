        C_b=FC{7};
        neigh_up=FC{12};
        neigh_down=FC{13};
        P=CELL{neigh_down(1,1)};
        Q=CELL{neigh_up(1,1)};
        C_j_down=norm_joint(C_b,FC{4}',P{5});
        C_j_up=norm_joint(C_b,FC{4}',Q{5});
       
        
        N=4900*1.5*2;
        tic;
        
        for k=1:N
            P1=[0;0];
            P2=[0;1];
            P3=[1;0];
            P=[0;5];
            x1=[0.5,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
            x2=[0.2,0,0,0,0,0,0,0,0];
            x3=[2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];
            x=triangle_2ndmapping(P1,P2,P3,x1,x2,x3,P);
        end
        toc;
        