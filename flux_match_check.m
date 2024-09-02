for q=1:9 % Change according to different pdf components
    FLX_ef=zeros(3,M); %%%% Edge flux
    FLX_efm=zeros(3,M); %%%% Edge flux match 1---match; 0---not match
    for r=1:M
        P=CELL{r};
        for i=1:3
            if i==1
                FLX_ef(i,r)=fl1(q,r);
            elseif i==2
                FLX_ef(i,r)=fl2(q,r);
            else
                FLX_ef(i,r)=fl3(q,r);
            end
        end
    end
    
    for r=1:M
        P=CELL{r};
        for i=1:3
            a=0; % The order number of the edge of the neighbor
            %%%% e.g. the ith edge of current triangle and the ath
            %%%% edge of its neighbor on that edge is the same edge
            if P{2+i-1}==0 %%%% No neighbor
                ;
            else
                Q=CELL{P{2+i-1}}; % The neighbor triangle
                for a=1:3
                    if Q{2+a-1}==r
                        break; %%%% Found a
                    end
                end
            end
            
            if a==0 %%%% No neighbor
                FLX_efm(i,r)=-1; %%%% No neighbor
            else
                if double(FLX_ef(i,r))==-double(FLX_ef(a,P{2+i-1}))
                    FLX_efm(i,r)=1; %%%% Flux matches
                else
                    ;
                end
            end
            
        end
    end
    
    for r=1:M
        for i=1:3
            if FLX_efm(i,r)==0
                error('dismatch');
            end
        end
    end
end