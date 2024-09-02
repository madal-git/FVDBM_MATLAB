P=0;
counter=0;
X=500;
for i=1:X
    for j=2:i-1
        if mod(i,j)==0
            break;
        else
            counter=counter+1;
            P(counter)=i;
        end
    end
end