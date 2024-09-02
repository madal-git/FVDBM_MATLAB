function f_2nd=in_cell_mapping_SLIC_apply(CELL,V,f,feq,f_nd,Tau,FPDC)

M=length(CELL);
[d,q]=size(V);
f_2nd=zeros(q,M);
feq_2nd=zeros(q,M);


%% Algotirhm 1, equally weighted 2nd-order interpolation See CELL v.11.11.2019
for i=1:M
    cell_center=CELL{i};
    %% Using three point, one of which is from boundary
    if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
        cell_neighbor=cell_center{57};
        CoeM=cell_center{58};
        F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
%         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
        f_2nd(:,i)=f(:,i)+diag(CoeM*F);
    elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
        error('Not available yet!');
%         cell_neighbor=cell_center{41};
%         fc_nd=cell_center{42};
%         CoeM=cell_center{44};
%         if cell_neighbor(1)==0
%             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%             f_1=2*fs_intcpt-f(:,s_cell);
% %             f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%             f_2=f(:,cell_neighbor(2));
%             f_3=f(:,cell_neighbor(3));
%         elseif cell_neighbor(2)==0
%             f_1=f(:,cell_neighbor(1));
%             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%             f_2=2*fs_intcpt-f(:,s_cell);
% %             f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%             f_3=f(:,cell_neighbor(3));
%         elseif cell_neighbor(3)==0
%             f_1=f(:,cell_neighbor(1));
%             f_2=f(:,cell_neighbor(2));
%             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%             f_3=2*fs_intcpt-f(:,s_cell);
% %             f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%         else
%             error('logic error!');
%         end
%         F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%         f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
% %         x=f(:,s_cell);
    elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
        if FPDC==0 % No periodic boundaries
            cell_neighbor=cell_center{59};
            fc_nd=cell_center{60};
            CoeM=cell_center{61};
            if cell_neighbor(1)==0
                f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_1=2*fs_intcpt-f(:,i);
                f_2=f(:,cell_neighbor(2));
                f_3=f(:,cell_neighbor(3));
            elseif cell_neighbor(2)==0
                f_1=f(:,cell_neighbor(1));
                f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_2=2*fs_intcpt-f(:,i);
                f_3=f(:,cell_neighbor(3));
            elseif cell_neighbor(3)==0
                f_1=f(:,cell_neighbor(1));
                f_2=f(:,cell_neighbor(2));
                f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_3=2*fs_intcpt-f(:,i);
            else
                error('logic error!');
            end
            F=[(f_1-f(:,i))';(f_2-f(:,i))';(f_3-f(:,i))'];
%             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
            f_2nd(:,i)=f(:,i)+diag(CoeM*F);
        elseif FPDC==1 % Only left & right boundaries are periodic
%             cell_neighbor=cell_center{49};
%             fc_nd=cell_center{50};
%             CoeM=cell_center{52};
%             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
%                 if cell_neighbor(1)==0
% %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_1=2*fs_intcpt-f(:,s_cell);
%                     f_2=f(:,cell_neighbor(2));
%                     f_3=f(:,cell_neighbor(3));
%                 elseif cell_neighbor(2)==0
%                     f_1=f(:,cell_neighbor(1));
% %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_2=2*fs_intcpt-f(:,s_cell);
%                     f_3=f(:,cell_neighbor(3));
%                 elseif cell_neighbor(3)==0
%                     f_1=f(:,cell_neighbor(1));
%                     f_2=f(:,cell_neighbor(2));
% %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_3=2*fs_intcpt-f(:,s_cell);
%                 else
%                     error('logic error!');
%                 end
%             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
        elseif FPDC==2 % Only top & bottom boundaries are periodic
%             cell_neighbor=cell_center{53};
%             fc_nd=cell_center{54};
%             CoeM=cell_center{56};
%             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
%                 if cell_neighbor(1)==0
% %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_1=2*fs_intcpt-f(:,s_cell);
%                     f_2=f(:,cell_neighbor(2));
%                     f_3=f(:,cell_neighbor(3));
%                 elseif cell_neighbor(2)==0
%                     f_1=f(:,cell_neighbor(1));
% %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_2=2*fs_intcpt-f(:,s_cell);
%                     f_3=f(:,cell_neighbor(3));
%                 elseif cell_neighbor(3)==0
%                     f_1=f(:,cell_neighbor(1));
%                     f_2=f(:,cell_neighbor(2));
% %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                     f_3=2*fs_intcpt-f(:,s_cell);
%                 else
%                     error('logic error!');
%                 end
% %                 x=f(:,s_cell);
%             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
%             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';   
        elseif FPDC==3 % All boundaries are periodic
            cell_neighbor=cell_center{62};
            CoeM=cell_center{64};
            F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
%             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
            f_2nd(:,i)=f(:,i)+diag(CoeM*F);
        else
            error('Logic error!');
        end
    else
        error('Logic error!');
    end
end



% %% Algotirhm 2, equally weighted 2nd-order interpolation + Relaxation
% for i=1:M
%     cell_center=CELL{i};
%     %% Using three point, one of which is from boundary
%     if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{57};
%         CoeM=cell_center{58};
%         F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
%         Feq=[(feq(:,cell_neighbor(1))-feq(:,i))';(feq(:,cell_neighbor(2))-feq(:,i))';(feq(:,cell_neighbor(3))-feq(:,i))'];
% %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%         f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         feq_2nd(:,i)=feq(:,i)+diag(CoeM*Feq);
%         f_2nd(:,i)=f_2nd(:,i)-(f_2nd(:,i)-feq_2nd(:,i))/Tau;
%     elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
% %         cell_neighbor=cell_center{41};
% %         fc_nd=cell_center{42};
% %         CoeM=cell_center{44};
% %         if cell_neighbor(1)==0
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_1=2*fs_intcpt-f(:,s_cell);
% % %             f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_2=f(:,cell_neighbor(2));
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(2)==0
% %             f_1=f(:,cell_neighbor(1));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_2=2*fs_intcpt-f(:,s_cell);
% % %             f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(3)==0
% %             f_1=f(:,cell_neighbor(1));
% %             f_2=f(:,cell_neighbor(2));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_3=2*fs_intcpt-f(:,s_cell);
% % %             f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %         else
% %             error('logic error!');
% %         end
% %         F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
% % %         x=f(:,s_cell);
%     elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
%         if FPDC==0 % No periodic boundaries
%             cell_neighbor=cell_center{59};
%             fc_nd=cell_center{60};
%             CoeM=cell_center{61};
%             if cell_neighbor(1)==0
% %                 f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_1=2*fs_intcpt-f(:,i);
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(2)==0
%                 f_1=f(:,cell_neighbor(1));
% %                 f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_2=2*fs_intcpt-f(:,i);
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(3)==0
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
% %                 f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_3=2*fs_intcpt-f(:,i);
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,i))';(f_2-f(:,i))';(f_3-f(:,i))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         elseif FPDC==1 % Only left & right boundaries are periodic
% %             cell_neighbor=cell_center{49};
% %             fc_nd=cell_center{50};
% %             CoeM=cell_center{52};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
%         elseif FPDC==2 % Only top & bottom boundaries are periodic
% %             cell_neighbor=cell_center{53};
% %             fc_nd=cell_center{54};
% %             CoeM=cell_center{56};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% % %                 x=f(:,s_cell);
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
%         elseif FPDC==3 % All boundaries are periodic
%             cell_neighbor=cell_center{62};
%             CoeM=cell_center{64};
%             F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
%             Feq=[(feq(:,cell_neighbor(1))-feq(:,i))';(feq(:,cell_neighbor(2))-feq(:,i))';(feq(:,cell_neighbor(3))-feq(:,i))'];
%             %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%             feq_2nd(:,i)=feq(:,i)+diag(CoeM*Feq);
%             f_2nd(:,i)=f_2nd(:,i)-(f_2nd(:,i)-feq_2nd(:,i))/Tau;
%         else
%             error('Logic error!');
%         end
%     else
%         error('Logic error!');
%     end
% end




% %% Algotirhm 3, upwind 2nd-order interpolation
% for i=1:M
%     cell_center=CELL{i};
%     %% Using three point, one of which is from boundary
%     if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{57};
%         CoeM=cell_center{58};
%         Zone_ID=cell_center{59};
% %         F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         for k=1:q
%            if Zone_ID(1,k)==123
%                f_2nd(k,i)=f(k,i);
%            elseif Zone_ID(1,k)==1
%                F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            elseif Zone_ID(1,k)==2
%                F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            elseif Zone_ID(1,k)==3
%                F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            elseif Zone_ID(1,k)==12
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%            elseif Zone_ID(1,k)==23
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%            elseif Zone_ID(1,k)==31
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%            else
%                error('Logic Error!');
%            end
%         end
%     elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
% %         cell_neighbor=cell_center{41};
% %         fc_nd=cell_center{42};
% %         CoeM=cell_center{44};
% %         if cell_neighbor(1)==0
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_1=2*fs_intcpt-f(:,s_cell);
% % %             f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_2=f(:,cell_neighbor(2));
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(2)==0
% %             f_1=f(:,cell_neighbor(1));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_2=2*fs_intcpt-f(:,s_cell);
% % %             f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(3)==0
% %             f_1=f(:,cell_neighbor(1));
% %             f_2=f(:,cell_neighbor(2));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_3=2*fs_intcpt-f(:,s_cell);
% % %             f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %         else
% %             error('logic error!');
% %         end
% %         F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
% % %         x=f(:,s_cell);
%     elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
%         if FPDC==0 % No periodic boundaries
%             cell_neighbor=cell_center{59};
%             fc_nd=cell_center{60};
%             CoeM=cell_center{61};
%             if cell_neighbor(1)==0
% %                 f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_1=2*fs_intcpt-f(:,i);
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(2)==0
%                 f_1=f(:,cell_neighbor(1));
% %                 f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_2=2*fs_intcpt-f(:,i);
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(3)==0
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
% %                 f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_3=2*fs_intcpt-f(:,i);
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,i))';(f_2-f(:,i))';(f_3-f(:,i))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         elseif FPDC==1 % Only left & right boundaries are periodic
% %             cell_neighbor=cell_center{49};
% %             fc_nd=cell_center{50};
% %             CoeM=cell_center{52};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
%         elseif FPDC==2 % Only top & bottom boundaries are periodic
% %             cell_neighbor=cell_center{53};
% %             fc_nd=cell_center{54};
% %             CoeM=cell_center{56};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% % %                 x=f(:,s_cell);
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';   
%         elseif FPDC==3 % All boundaries are periodic
%             cell_neighbor=cell_center{64};
%             CoeM=cell_center{66};
%             Zone_ID=cell_center{67};
%             
%             for k=1:q
%                 if Zone_ID(1,k)==123
%                     f_2nd(k,i)=f(k,i);
%                 elseif Zone_ID(1,k)==1
%                     F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 elseif Zone_ID(1,k)==2
%                     F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 elseif Zone_ID(1,k)==3
%                     F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 elseif Zone_ID(1,k)==12
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%                 elseif Zone_ID(1,k)==23
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%                 elseif Zone_ID(1,k)==31
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%                 else
%                     error('Logic Error!');
%                 end
%             end
%         else
%             error('Logic error!');
%         end
%     else
%         error('Logic error!');
%     end
% end


% %% Algotirhm 4, upwind 2nd-order interpolation + Relaxation at the SLIC points
% for i=1:M
%     cell_center=CELL{i};
%     %% Using three point, one of which is from boundary
%     if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{57};
%         CoeM=cell_center{58};
%         Zone_ID=cell_center{59};
% %         F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         for k=1:q
%            if Zone_ID(1,k)==123
%                f_2nd(k,i)=f(k,i);
%                feq_2nd(k,i)=feq(k,i);
%            elseif Zone_ID(1,k)==1
%                F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                Feq=[(feq(k,cell_neighbor(1))-feq(k,i)),(feq(k,cell_neighbor(2))-feq(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%            elseif Zone_ID(1,k)==2
%                F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                Feq=[(feq(k,cell_neighbor(2))-feq(k,i)),(feq(k,cell_neighbor(3))-feq(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%            elseif Zone_ID(1,k)==3
%                F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                Feq=[(feq(k,cell_neighbor(3))-feq(k,i)),(feq(k,cell_neighbor(1))-feq(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%            elseif Zone_ID(1,k)==12
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%                feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(2))];
%            elseif Zone_ID(1,k)==23
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%                feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(3))];
%            elseif Zone_ID(1,k)==31
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%                feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(1))];
%            else
%                error('Logic Error!');
%            end
%         end
%         f_2nd(:,i)=f_2nd(:,i)-(f_2nd(:,i)-feq_2nd(:,i))/Tau;
%     elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
% %         cell_neighbor=cell_center{41};
% %         fc_nd=cell_center{42};
% %         CoeM=cell_center{44};
% %         if cell_neighbor(1)==0
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_1=2*fs_intcpt-f(:,s_cell);
% % %             f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_2=f(:,cell_neighbor(2));
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(2)==0
% %             f_1=f(:,cell_neighbor(1));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_2=2*fs_intcpt-f(:,s_cell);
% % %             f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(3)==0
% %             f_1=f(:,cell_neighbor(1));
% %             f_2=f(:,cell_neighbor(2));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_3=2*fs_intcpt-f(:,s_cell);
% % %             f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %         else
% %             error('logic error!');
% %         end
% %         F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
% % %         x=f(:,s_cell);
%     elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
%         if FPDC==0 % No periodic boundaries
%             cell_neighbor=cell_center{59};
%             fc_nd=cell_center{60};
%             CoeM=cell_center{61};
%             if cell_neighbor(1)==0
% %                 f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_1=2*fs_intcpt-f(:,i);
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(2)==0
%                 f_1=f(:,cell_neighbor(1));
% %                 f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_2=2*fs_intcpt-f(:,i);
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(3)==0
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
% %                 f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_3=2*fs_intcpt-f(:,i);
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,i))';(f_2-f(:,i))';(f_3-f(:,i))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         elseif FPDC==1 % Only left & right boundaries are periodic
% %             cell_neighbor=cell_center{49};
% %             fc_nd=cell_center{50};
% %             CoeM=cell_center{52};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
%         elseif FPDC==2 % Only top & bottom boundaries are periodic
% %             cell_neighbor=cell_center{53};
% %             fc_nd=cell_center{54};
% %             CoeM=cell_center{56};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% % %                 x=f(:,s_cell);
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';   
%         elseif FPDC==3 % All boundaries are periodic
%             cell_neighbor=cell_center{64};
%             CoeM=cell_center{66};
%             Zone_ID=cell_center{67};
%             
%             for k=1:q
%                 if Zone_ID(1,k)==123
%                     f_2nd(k,i)=f(k,i);
%                     feq_2nd(k,i)=feq(k,i);
%                 elseif Zone_ID(1,k)==1
%                     F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                     Feq=[(feq(k,cell_neighbor(1))-feq(k,i)),(feq(k,cell_neighbor(2))-feq(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                     feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%                 elseif Zone_ID(1,k)==2
%                     F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                     Feq=[(feq(k,cell_neighbor(2))-feq(k,i)),(feq(k,cell_neighbor(3))-feq(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                     feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%                 elseif Zone_ID(1,k)==3
%                     F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                     Feq=[(feq(k,cell_neighbor(3))-feq(k,i)),(feq(k,cell_neighbor(1))-feq(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                     feq_2nd(k,i)=feq(k,i)+(CoeM(k,:)*Feq)';
%                 elseif Zone_ID(1,k)==12
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%                     feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(2))];
%                 elseif Zone_ID(1,k)==23
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%                     feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(3))];
%                 elseif Zone_ID(1,k)==31
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%                     feq_2nd(k,i)=CoeM(k,:)*[feq(k,i);feq(k,cell_neighbor(1))];
%                 else
%                     error('Logic Error!');
%                 end
%             end
%             f_2nd(:,i)=f_2nd(:,i)-(f_2nd(:,i)-feq_2nd(:,i))/Tau;
%         else
%             error('Logic error!');
%         end
%     else
%         error('Logic error!');
%     end
% end


%% Algotirhm 5, downwind 2nd-order interpolation
% for i=1:M
%     cell_center=CELL{i};
%     %% Using three point, one of which is from boundary
%     if cell_center{37}==1 % Type-1 cell, see CELL v.06.16.2017
%         cell_neighbor=cell_center{57};
%         CoeM=cell_center{58};
%         Zone_ID=cell_center{59};
% %         F=[(f(:,cell_neighbor(1))-f(:,i))';(f(:,cell_neighbor(2))-f(:,i))';(f(:,cell_neighbor(3))-f(:,i))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd(:,i)=f(:,i)+diag(CoeM*F); 
%         for k=1:q
%            if Zone_ID(1,k)==123
%                f_2nd(k,i)=f(k,i);
%            elseif Zone_ID(1,k)==1
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%            elseif Zone_ID(1,k)==2
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%            elseif Zone_ID(1,k)==3
%                f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%            elseif Zone_ID(1,k)==12
%                F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            elseif Zone_ID(1,k)==23
%                F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            elseif Zone_ID(1,k)==31
%                F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%            else
%                error('Logic Error!');
%            end
%         end
%     elseif cell_center{37}==3 % Type-3 cell, see CELL v.06.16.2017
% %         cell_neighbor=cell_center{41};
% %         fc_nd=cell_center{42};
% %         CoeM=cell_center{44};
% %         if cell_neighbor(1)==0
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_1=2*fs_intcpt-f(:,s_cell);
% % %             f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_2=f(:,cell_neighbor(2));
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(2)==0
% %             f_1=f(:,cell_neighbor(1));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_2=2*fs_intcpt-f(:,s_cell);
% % %             f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %             f_3=f(:,cell_neighbor(3));
% %         elseif cell_neighbor(3)==0
% %             f_1=f(:,cell_neighbor(1));
% %             f_2=f(:,cell_neighbor(2));
% %             fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %             f_3=2*fs_intcpt-f(:,s_cell);
% % %             f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %         else
% %             error('logic error!');
% %         end
% %         F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %         x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %         f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
% % %         x=f(:,s_cell);
%     elseif cell_center{37}==2 % Type-2 cell, see CELL v.06.16.2017
%         if FPDC==0 % No periodic boundaries
%             cell_neighbor=cell_center{59};
%             fc_nd=cell_center{60};
%             CoeM=cell_center{61};
%             if cell_neighbor(1)==0
% %                 f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_1=2*fs_intcpt-f(:,i);
%                 f_2=f(:,cell_neighbor(2));
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(2)==0
%                 f_1=f(:,cell_neighbor(1));
% %                 f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_2=2*fs_intcpt-f(:,i);
%                 f_3=f(:,cell_neighbor(3));
%             elseif cell_neighbor(3)==0
%                 f_1=f(:,cell_neighbor(1));
%                 f_2=f(:,cell_neighbor(2));
% %                 f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
%                 fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
%                 f_3=2*fs_intcpt-f(:,i);
%             else
%                 error('logic error!');
%             end
%             F=[(f_1-f(:,i))';(f_2-f(:,i))';(f_3-f(:,i))'];
% %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
%             f_2nd(:,i)=f(:,i)+diag(CoeM*F);
%         elseif FPDC==1 % Only left & right boundaries are periodic
% %             cell_neighbor=cell_center{49};
% %             fc_nd=cell_center{50};
% %             CoeM=cell_center{52};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';
%         elseif FPDC==2 % Only top & bottom boundaries are periodic
% %             cell_neighbor=cell_center{53};
% %             fc_nd=cell_center{54};
% %             CoeM=cell_center{56};
% %             if length(setxor(cell_neighbor,0))==2 % The current cell is attached to a non-periodic boundary
% %                 if cell_neighbor(1)==0
% % %                     f_1=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_1=2*fs_intcpt-f(:,s_cell);
% %                     f_2=f(:,cell_neighbor(2));
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(2)==0
% %                     f_1=f(:,cell_neighbor(1));
% % %                     f_2=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_2=2*fs_intcpt-f(:,s_cell);
% %                     f_3=f(:,cell_neighbor(3));
% %                 elseif cell_neighbor(3)==0
% %                     f_1=f(:,cell_neighbor(1));
% %                     f_2=f(:,cell_neighbor(2));
% % %                     f_3=(f_nd(:,fc_nd(1))+f_nd(:,fc_nd(2)))/2;
% %                     fs_intcpt=(1-fc_nd(3))*f_nd(:,fc_nd(1))+fc_nd(3)*f_nd(:,fc_nd(2));
% %                     f_3=2*fs_intcpt-f(:,s_cell);
% %                 else
% %                     error('logic error!');
% %                 end
% % %                 x=f(:,s_cell);
% %             elseif length(setxor(cell_neighbor,0))==4 % The current cell is attached to a periodic boundary
% %                 f_1=f(:,cell_neighbor(1));
% %                 f_2=f(:,cell_neighbor(2));
% %                 f_3=f(:,cell_neighbor(3));
% %             else
% %                 error('logic error!');
% %             end
% %             F=[(f_1-f(:,s_cell))';(f_2-f(:,s_cell))';(f_3-f(:,s_cell))'];
% % %             x=f(:,s_cell)+((s_stcl_coord-cell_center{5})'/sqrt(cell_center{6})*(CoeA*F))';      
% %             f_2nd=f(:,s_cell)+((s_stcl_coord-cell_center{5})'*(CoeM*F))';   
%         elseif FPDC==3 % All boundaries are periodic
%             cell_neighbor=cell_center{64};
%             CoeM=cell_center{66};
%             Zone_ID=cell_center{67};
%             
%             for k=1:q
%                 if Zone_ID(1,k)==123
%                     f_2nd(k,i)=f(k,i);
%                 elseif Zone_ID(1,k)==1
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(3))];
%                 elseif Zone_ID(1,k)==2
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(1))];
%                 elseif Zone_ID(1,k)==3
%                     f_2nd(k,i)=CoeM(k,:)*[f(k,i);f(k,cell_neighbor(2))];
%                 elseif Zone_ID(1,k)==12
%                     F=[(f(k,cell_neighbor(3))-f(k,i)),(f(k,cell_neighbor(1))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 elseif Zone_ID(1,k)==23
%                     F=[(f(k,cell_neighbor(1))-f(k,i)),(f(k,cell_neighbor(2))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 elseif Zone_ID(1,k)==31
%                     F=[(f(k,cell_neighbor(2))-f(k,i)),(f(k,cell_neighbor(3))-f(k,i))]';
%                     f_2nd(k,i)=f(k,i)+(CoeM(k,:)*F)';
%                 else
%                     error('Logic Error!');
%                 end
%             end
%         else
%             error('Logic error!');
%         end
%     else
%         error('Logic error!');
%     end
% end