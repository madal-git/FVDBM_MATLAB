function [x_2nd,Cxy]=pfls(coord,current_cell_num,CELL,x_cl,x_nd,FPDC)

% [x_2nd,Cxy]=pfls(coord,current_cell_num,CELL,x_cl,x_nd,FPDC) calculates the
% value at any location with the plane-fitting least square (PFLS) method,
% which is 2nd-order in space as well as the c1, c2 coefficients in 2nd-order reconstruction of f = c1*x+c2*y+c3
cell_center=CELL{current_cell_num};
%% Using three point, one of which is from boundary
if cell_center{37}==1 % Type-1 cell
    cell_neighbor=cell_center{38};
    CoeA=cell_center{40};
    F=[(x_cl(:,cell_neighbor(1))-x_cl(:,current_cell_num)),(x_cl(:,cell_neighbor(2))-x_cl(:,current_cell_num)),(x_cl(:,cell_neighbor(3))-x_cl(:,current_cell_num))];
    Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
    x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
elseif cell_center{37}==3 % Type-3 cell
    cell_neighbor=cell_center{41};
    fc_nd=cell_center{42};
    CoeA=cell_center{44};
    if cell_neighbor(1)==0
        fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
        x_c_1=2*fs_intcpt-x_cl(:,current_cell_num);
        %             x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
        x_c_2=x_cl(:,cell_neighbor(2));
        x_c_3=x_cl(:,cell_neighbor(3));
    elseif cell_neighbor(2)==0
        x_c_1=x_cl(:,cell_neighbor(1));
        fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
        x_c_2=2*fs_intcpt-x_cl(:,current_cell_num);
        %             x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
        x_c_3=x_cl(:,cell_neighbor(3));
    elseif cell_neighbor(3)==0
        x_c_1=x_cl(:,cell_neighbor(1));
        x_c_2=x_cl(:,cell_neighbor(2));
        fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
        x_c_3=2*fs_intcpt-x_cl(:,current_cell_num);
        %             x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
    else
        error('logic error!');
    end
    F=[(x_c_1-x_cl(:,current_cell_num)),(x_c_2-x_cl(:,current_cell_num)),(x_c_3-x_cl(:,current_cell_num))];
    Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
    x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
elseif cell_center{37}==2 % Type-2 cell
    if FPDC==0 % No periodic boundaries
        cell_neighbor=cell_center{41};
        fc_nd=cell_center{42};
        CoeA=cell_center{44};
        if cell_neighbor(1)==0
            %                 x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_1=2*fs_intcpt-x_cl(:,current_cell_num);
            x_c_2=x_cl(:,cell_neighbor(2));
            x_c_3=x_cl(:,cell_neighbor(3));
        elseif cell_neighbor(2)==0
            x_c_1=x_cl(:,cell_neighbor(1));
            %                 x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_2=2*fs_intcpt-x_cl(:,current_cell_num);
            x_c_3=x_cl(:,cell_neighbor(3));
        elseif cell_neighbor(3)==0
            x_c_1=x_cl(:,cell_neighbor(1));
            x_c_2=x_cl(:,cell_neighbor(2));
            %                 x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
            fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
            x_c_3=2*fs_intcpt-x_cl(:,current_cell_num);
        else
            error('logic error!');
        end
        F=[(x_c_1-x_cl(:,current_cell_num)),(x_c_2-x_cl(:,current_cell_num)),(x_c_3-x_cl(:,current_cell_num))];
        Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
        x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
    elseif FPDC==1 % Only left & right boundaries are periodic
        cell_neighbor=cell_center{49};
        fc_nd=cell_center{50};
        CoeA=cell_center{52};
        if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
            if cell_neighbor(1)==0
                %                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_1=2*fs_intcpt-x_cl(:,current_cell_num);
                x_c_2=x_cl(:,cell_neighbor(2));
                x_c_3=x_cl(:,cell_neighbor(3));
            elseif cell_neighbor(2)==0
                x_c_1=x_cl(:,cell_neighbor(1));
                %                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_2=2*fs_intcpt-x_cl(:,current_cell_num);
                x_c_3=x_cl(:,cell_neighbor(3));
            elseif cell_neighbor(3)==0
                x_c_1=x_cl(:,cell_neighbor(1));
                x_c_2=x_cl(:,cell_neighbor(2));
                %                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_3=2*fs_intcpt-x_cl(:,current_cell_num);
            else
                error('logic error!');
            end
        elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
            x_c_1=x_cl(:,cell_neighbor(1));
            x_c_2=x_cl(:,cell_neighbor(2));
            x_c_3=x_cl(:,cell_neighbor(3));
        else
            error('logic error!');
        end
        F=[(x_c_1-x_cl(:,current_cell_num)),(x_c_2-x_cl(:,current_cell_num)),(x_c_3-x_cl(:,current_cell_num))];
        Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
        x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
    elseif FPDC==2 % Only top & bottom boundaries are periodic
        cell_neighbor=cell_center{53};
        fc_nd=cell_center{54};
        CoeA=cell_center{56};
        if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
            if cell_neighbor(1)==0
                %                     x_c_1=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_1=2*fs_intcpt-x_cl(:,current_cell_num);
                x_c_2=x_cl(:,cell_neighbor(2));
                x_c_3=x_cl(:,cell_neighbor(3));
            elseif cell_neighbor(2)==0
                x_c_1=x_cl(:,cell_neighbor(1));
                %                     x_c_2=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_2=2*fs_intcpt-x_cl(:,current_cell_num);
                x_c_3=x_cl(:,cell_neighbor(3));
            elseif cell_neighbor(3)==0
                x_c_1=x_cl(:,cell_neighbor(1));
                x_c_2=x_cl(:,cell_neighbor(2));
                %                     x_c_3=(x_nd(:,fc_nd(1))+x_nd(:,fc_nd(2)))/2;
                fs_intcpt=(1-fc_nd(3))*x_nd(:,fc_nd(1))+fc_nd(3)*x_nd(:,fc_nd(2));
                x_c_3=2*fs_intcpt-x_cl(:,current_cell_num);
            else
                error('logic error!');
            end
            %                 x_2nd=x_cl(:,current_cell_num);
        elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
            x_c_1=x_cl(:,cell_neighbor(1));
            x_c_2=x_cl(:,cell_neighbor(2));
            x_c_3=x_cl(:,cell_neighbor(3));
        else
            error('logic error!');
        end
        F=[(x_c_1-x_cl(:,current_cell_num)),(x_c_2-x_cl(:,current_cell_num)),(x_c_3-x_cl(:,current_cell_num))];
        Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
        x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
    elseif FPDC==3 % All boundaries are periodic
        cell_neighbor=cell_center{45};
        CoeA=cell_center{48};
        F=[(x_cl(:,cell_neighbor(1))-x_cl(:,current_cell_num)),(x_cl(:,cell_neighbor(2))-x_cl(:,current_cell_num)),(x_cl(:,cell_neighbor(3))-x_cl(:,current_cell_num))];
        Cxy=F*CoeA; % This is the c1 and c2 coefficients in f = c1*x+c2*y+c3, which is also the df/dx and df/dy that will be used to compute the Oscillation Indicator
        x_2nd=x_cl(:,current_cell_num)+Cxy*(coord-cell_center{5});
    else
        error('Logic error!');
    end
else
    error('Logic error!');
end