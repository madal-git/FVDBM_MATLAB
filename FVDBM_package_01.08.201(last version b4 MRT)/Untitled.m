ST=0;
for l=1:O
    FC=FACE{l};
    St=FC{16};
    S1=St{10};
    S2=St{11};
    S=unique(union(S1,S2));
    if length(S)>2
        error('?');
    end
    ST=union(ST,S);
end
if length(ST)-1~=length(setxor(ST,0))
    error('!');
end
ST=setxor(ST,0);
unique(ST)

for l=1:O
    FC=FACE{l};
    St=FC{16};
    Cd=St{1};
    Cu=St{2};
    Sd=St{10};
    Su=St{11};
    xd=St{4};
    yd=St{5};
    xu=St{6};
    yu=St{7};
    L1=length(Cd);
    L2=length(Cu);
    L3=length(Sd);
    L4=length(Su);
    L5=length(xd);
    L6=length(yd);
    L7=length(xu);
    L8=length(yu);
    if (L1~=qh || L2~=qh || L3~=qh || L4~=qh) || (L5~=qh || L6~=qh || L7~=qh || L8~=qh)
        error('h');
    end
    if FC{2}==0
        for k=1:qh
            Cell_down=CELL{Cd(1,k)};
            Cell_up=CELL{Cu(1,k)};
            nd_down_1=Cell_down{13};
            nd_down_2=Cell_down{14};
            nd_down_3=Cell_down{15};
            nd_up_1=Cell_up{13};
            nd_up_2=Cell_up{14};
            nd_up_3=Cell_up{15};
            in_zone_one=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_1,nd_down_2);
            in_zone_two=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_2,nd_down_3);
            in_zone_three=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_3,nd_down_1);
            if in_zone_one
                if Sd(1,k)~=1
                    error('1');
                end% Located in zone 1
            elseif in_zone_two
                if Sd(1,k)~=2
                    error('2');
                end% 
            elseif in_zone_three
                if Sd(1,k)~=3
                    error('3');
                end% 
            else
                on_edge_one=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_1);
                on_edge_two=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_2);
                on_edge_three=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_3);
                if on_edge_one && on_edge_two && on_edge_three
                    if Sd(1,k)~=4
                        error('4');
                    end%
                elseif on_edge_one
                    if Sd(1,k)~=1
                        error('1');
                    end%
                elseif on_edge_two
                    if Sd(1,k)~=2
                        error('2');
                    end%
                elseif on_edge_three
                    if Sd(1,k)~=3
                        error('3');
                    end%
                else
                    error('The zone for the stencil point is not found!');
                end
            end
            % upwind
            in_zone_one=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_1,nd_up_2);
            in_zone_two=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_2,nd_up_3);
            in_zone_three=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_3,nd_up_1);
            if in_zone_one
                if Su(1,k)~=1
                    error('1');
                end% Located in zone 1
            elseif in_zone_two
                if Su(1,k)~=2
                    error('2');
                end% 
            elseif in_zone_three
                if Su(1,k)~=3
                    error('3');
                end% 
            else
                on_edge_one=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_1);
                on_edge_two=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_2);
                on_edge_three=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_3);
                if on_edge_one && on_edge_two && on_edge_three
                    if Su(1,k)~=4
                        error('4');
                    end%
                elseif on_edge_one
                    if Su(1,k)~=1
                        error('1');
                    end%
                elseif on_edge_two
                    if Su(1,k)~=2
                        error('2');
                    end%
                elseif on_edge_three
                    if Su(1,k)~=3
                        error('3');
                    end%
                else
                    error('The zone for the stencil point is not found!');
                end
            end
        end
    else
        for k=1:qh
            if Cd(1,k)==0 && Cu(1,k)~=0
                Cell_up=CELL{Cu(1,k)};
                nd_up_1=Cell_up{13};
                nd_up_2=Cell_up{14};
                nd_up_3=Cell_up{15};
                if Sd(1,k)~=4
                    error('4');
                end
                in_zone_one=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_1,nd_up_2);
                in_zone_two=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_2,nd_up_3);
                in_zone_three=in_triangle([xu(1,k);yu(1,k)],Cell_up{5},nd_up_3,nd_up_1);
                if in_zone_one
                    if Su(1,k)~=1
                        error('1');
                    end% Located in zone 1
                elseif in_zone_two
                    if Su(1,k)~=2
                        error('2');
                    end%
                elseif in_zone_three
                    if Su(1,k)~=3
                        error('3');
                    end%
                else
                    on_edge_one=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_1);
                    on_edge_two=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_2);
                    on_edge_three=on_edge([xu(1,k);yu(1,k)],Cell_up{5},nd_up_3);
                    if on_edge_one && on_edge_two && on_edge_three
                        if Su(1,k)~=4
                            error('4');
                        end%
                    elseif on_edge_one
                        if Su(1,k)~=1
                            error('1');
                        end%
                    elseif on_edge_two
                        if Su(1,k)~=2
                            error('2');
                        end%
                    elseif on_edge_three
                        if Su(1,k)~=3
                            error('3');
                        end%
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
            elseif Cd(1,k)~=0 && Cu(1,k)==0
                Cell_down=CELL{Cd(1,k)};
                nd_down_1=Cell_down{13};
                nd_down_2=Cell_down{14};
                nd_down_3=Cell_down{15};
                if Su(1,k)~=4
                    error('4');
                end
                in_zone_one=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_1,nd_down_2);
                in_zone_two=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_2,nd_down_3);
                in_zone_three=in_triangle([xd(1,k);yd(1,k)],Cell_down{5},nd_down_3,nd_down_1);
                if in_zone_one
                    if Sd(1,k)~=1
                        error('1');
                    end% Located in zone 1
                elseif in_zone_two
                    if Sd(1,k)~=2
                        error('2');
                    end%
                elseif in_zone_three
                    if Sd(1,k)~=3
                        error('3');
                    end%
                else
                    on_edge_one=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_1);
                    on_edge_two=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_2);
                    on_edge_three=on_edge([xd(1,k);yd(1,k)],Cell_down{5},nd_down_3);
                    if on_edge_one && on_edge_two && on_edge_three
                        if Sd(1,k)~=4
                            error('4');
                        end%
                    elseif on_edge_one
                        if Sd(1,k)~=1
                            error('1');
                        end%
                    elseif on_edge_two
                        if Sd(1,k)~=2
                            error('2');
                        end%
                    elseif on_edge_three
                        if Sd(1,k)~=3
                            error('3');
                        end%
                    else
                        error('The zone for the stencil point is not found!');
                    end
                end
            else
                error('L');
            end
        end
    end
end