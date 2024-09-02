function v_n=norm_v(v)
% function v_n=norm_v(v) change the length of a 1D, 2D or 3D vector/vectors to unit while keep
% its/their direction unchanged


[D,L]=size(v);

if D>3 || D<=0
    error('The vector(s) must be 1D, 2D or 3D!');
end

if L==0
    error('There is no vector to be normalized!');
end

v_n=v;

for i=1:L
    if D==1
        l=v_n(1,i);
    elseif D==2
        l=sqrt(v_n(1,i)^2+v_n(2,i)^2);
    elseif D==3
        l=sqrt(v_n(1,i)^2+v_n(2,i)^2+v_n(3,i)^2);
    else
        error('The vector(s) must be 1D, 2D or 3D!');
    end
    if single(l+10)==single(10)
        ;
    else
        v_n(:,i)=v_n(:,i)/l;
    end
end