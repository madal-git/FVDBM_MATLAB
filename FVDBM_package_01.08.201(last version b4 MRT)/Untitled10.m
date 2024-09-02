        FC=FACE{l};
        fc23=FC{23};
        if fc23==2 || fc23==3 % Face type 2 and 3
            disp('1')
            fc16=FC{16};
            fc24=FC{24};
            fc26=FC{26};
            fc27=FC{20};
            for i=1:length(fc16)
                disp('2')
                if fc16(1,i)==0
                    disp('3')
                    if fc26(1,i)==4
                        fc27(1,i)=0;
                        fc20(2,i)=fc24(1,i);
                        fc20(3,i)=fc24(1,i);
                        fc20(4,i)=fc24(1,i);
                    else
                        disp('4')
                        Cell_fix=CELL{fc24(1,i)};
                        Face_target=FACE{Cell_fix{15+fc26(1,i)}};
                        Face_target{10}
                        Face_target{11}
                        if (Face_target{10}~=0 && Face_target{11}~=0)% || (Face_target{10}==0 && Face_target{11}~=0) % Both end nodes of target face are located on boundary
                            disp('5')
                            error('asdfasdf');
                        end
                    end
                end
            end
        end
        
                                    cell_pool_2_2_2=[10 15];
                            Cell_pool_2_2_2=cell(2,length(cell_pool_2_2_2));
                            for i=1:length(cell_pool_2_2_2)
                                C=Cell_pool_2_2_2{i};
                                Cell_p=CELL{cell_pool_2_2_2(i)};
                                C{1,1}=Cell_p{1};
                                C{2,1}=Cell_p{5}+[-(X2-X1);0];
                                Cell_pool_2_2_2{i}=C;
                            end