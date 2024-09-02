classdef tracer < handle
    properties
        L_x;
        L_y;
        num_x;
        num_y;
        N;
        location_old;
        location_new;
        velocity_old;
        velocity_new;
        current_cell_num_old;
        current_cell_num_new;
        FIG;
        FIG_COUNTER;
        e=1e-3;
        %% External variables channeled in
        X1;
        X2;
        Y1;
        Y2;
        U;
        U_nd;
        CELL;
        M;
        NODE;
        dt;
        FPDC;
    end
    
    methods
        function ini(obj)
            dx=obj.L_x/(obj.num_x+1);
            dy=obj.L_y/(obj.num_y+1);
            obj.location_old={[0;0]};
            obj.velocity_old={[0;0]};
            obj.N=obj.num_x*obj.num_y;
            obj.location_old=cell(1,obj.N);
            obj.location_new=cell(1,obj.N);
            obj.velocity_old=cell(1,obj.N);
            obj.velocity_new=cell(1,obj.N);
            obj.current_cell_num_old=cell(1,obj.N);
            obj.current_cell_num_new=cell(1,obj.N);
            counter=0;
            for i=1:obj.num_x
                for j=1:obj.num_y
                    counter=counter+1;
                    obj.location_old{counter}=[dx*i;dy*j];
                    obj.location_new{counter}=[dx*i;dy*j];
                end
            end
            if counter~=obj.N
                error('Logic Error!');
            end
            for i=1:obj.N
                obj.velocity_old{i}=[0;0];
                obj.velocity_new{i}=[0;0];
                obj.current_cell_num_old{i}=0;
                obj.current_cell_num_new{i}=0;
            end
            
            % Find the cells that contains each tracer point
            for i=1:obj.N
                r=in_which_triangle(obj.location_old{i},obj.CELL,obj.M);
                if r==0 % The point is not within any triangle
                    [nd1,nd2,ratio]=on_which_edge(obj.location_old{i},obj.CELL,obj.M);
                    if nd1==0 || nd2==0 % The given point is not on any edge
                        error('Logic Error!');
                    else
                        ND1=obj.NODE{nd1};
                        ND2=obj.NODE{nd2};
                        cell_common=intersect(ND1{5},ND2{5});
                        if cell_common==0
                            error('Logic Error!');
                        else
                            obj.current_cell_num_old{i}=cell_common(1);
                            obj.current_cell_num_new{i}=obj.current_cell_num_old{i};
                        end
                    end
                else
                    obj.current_cell_num_old{i}=r;
                    obj.current_cell_num_new{i}=obj.current_cell_num_old{i};
                end
            end
            
            mkdir tracerplot1
            obj.FIG_COUNTER=0;
            pathname = 'G:\Dropbox\Pitt!\LBM\repos\fvdbm\tracerplot1'; % Change this for different platform
            filename =['TracerPlot',num2str(obj.FIG_COUNTER,'%04d')];
            obj.show;
            saveas(gca,fullfile(pathname, filename),'jpg')
            print('-dtiff','-r300',filename);
        end
        
        function update(obj)
            for i=1:obj.N
                %% New velocity
                [obj.velocity_new{i},obj.current_cell_num_new{i}]=point_value(obj.location_old{i},obj.current_cell_num_old{i},obj.U,obj.U_nd,obj.CELL,obj.M,obj.NODE,obj.FPDC);
                %% New location
                obj.location_new{i}=obj.location_old{i}+(obj.velocity_old{i}+obj.velocity_new{i})/2*obj.dt;
                %% Avoid stuck at the top-right corner for lid-driven cavity flow
                Coord=obj.location_new{i};
                if dis(Coord,[obj.X2;obj.Y2])<=2*obj.e
                    Coord=[obj.X2;obj.Y2]-[1;1]*obj.e*5;
%                     obj.velocity_new{i}=norm(obj.velocity_new{i})*100*[-1;-1]; % bounce velocity
                    DIS=zeros(1,obj.M);
                    for j=1:obj.M
                        CL=obj.CELL{j};
                        DIS(j)=dis(CL{5},Coord);
                    end
                    DIS_min=min(DIS);
                    for j=1:obj.M
                        if single(DIS(j))==single(DIS_min)
                            break;
                        end
                    end
                    if j==obj.M
                        if DIS(i)~=DIS_min
                            error('Logic Error!');
                        end
                    end
                    obj.velocity_new{i}=obj.U(:,j);
                end
                %% Avoid stuck on the boundaries
                if Coord(1,1)<=obj.X1+obj.e
                    Coord(1,1)=obj.X1+obj.e;
                elseif Coord(1,1)>=obj.X2-obj.e
                    Coord(1,1)=obj.X2-obj.e;
                end
                if Coord(2,1)<=obj.Y1+obj.e
                    Coord(2,1)=obj.Y1+obj.e;
                elseif Coord(2,1)>=obj.Y2-obj.e
                    Coord(2,1)=obj.Y2-obj.e;
                end
                obj.location_new{i}=Coord;
                %% Final update
                obj.location_old{i}=obj.location_new{i};
                obj.velocity_old{i}=obj.velocity_new{i};
                obj.current_cell_num_old{i}=obj.current_cell_num_new{i};
            end
        end
        
        function show(obj)
            obj.FIG=figure(135);
            clf;
            plt_line([0;0],[obj.L_x;0]);
            hold on;
            plt_line([obj.L_x;0],[obj.L_x;obj.L_y]);
            hold on;
            plt_line([obj.L_x;obj.L_y],[0;obj.L_y]);
            hold on;
            plt_line([0;obj.L_y],[0;0]);
            hold on;
            for i=1:obj.N
                plt_point(obj.location_new{i});
                hold on;
            end
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            set(gcf,'color','w');
        end
        
        function showsave(obj)
            pathname = 'G:\Dropbox\Pitt!\LBM\repos\fvdbm\tracerplot1'; % Change this for different platform
            filename =['TracerPlot',num2str(obj.FIG_COUNTER,'%04d')];
            obj.show;
            saveas(gca,fullfile(pathname, filename),'jpg')
            print('-dtiff','-r300',filename)
        end
    end
end