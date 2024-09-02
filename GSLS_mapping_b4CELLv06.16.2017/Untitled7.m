for l=1:O
    FC=FACE{l};
    if FC{10}==0 && FC{11}==0
        if norm(FC{20})==0
            error(':');
        end
    end
end