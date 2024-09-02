function [nd1,nd2,ratio]=on_which_edge(p,CELL,M)
% function [nd1,nd2]=on_whcih_edge(p,CELL,M,ref) find on which edge the given
% point is located among all the edges of all triangles
% p is the coordinates of given point, has to be column vector
% CELL is the cell data structure of all triangles
% M is the total number of triangles
% ref is the reference length for killing roundoff error
% nd1 is the number of first node of the edge
% nd2 is the number of the second node of the edge
% ratio is the the length ratio intersected by the given point, ratio =
% dis(p,nd1)/dis(nd1,nd2)

a=0;
for r=1:M
    P=CELL{r};
    for s=1:3
        if s==1
            Ndl=P{13};
            Ndr=P{14};
        elseif s==2
            Ndl=P{14};
            Ndr=P{15};
        else
            Ndl=P{15};
            Ndr=P{13};
        end
        if on_edge(p,Ndl,Ndr)
            a=1;
            if s==1
                nd1=P{7};
                nd2=P{8};
            elseif s==2
                nd1=P{8};
                nd2=P{9};
            else
                nd1=P{9};
                nd2=P{7};
            end
            ratio=dis(p,Ndl)/dis(Ndl,Ndr);
            break;
        end
    end
    if a==1
        break;
    end
end

if a==0 % The given point is not on any edge
    nd1=0;
    nd2=0;
    ratio=0;
else
    if ratio==0 % The given point is located exactly at the left end point of found edge
        nd2=nd1;
    elseif ratio==1 % The given point is located exactly at the right end point of found edge
        nd1=nd2;
    else % % The given point is located in the middle of found edge
        ;
    end
end