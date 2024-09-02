function ang=angle(u,v)
% function ang=angle(u,v) determines the angle in degree between two vectors
% u and v in the x-y plane.
% where u=u1*i + u2*j, v=v1*i+v2*j and has to be column vector
% ang is the angle between them, which is always smaller than 90 degrees
e=1e-8;
u=u/norm(u);
v=v/norm(v);
L=abs(u(1,1)*v(2,1)-u(2,1)*v(1,1));
if L>1
    if (L-1)<e
        L=1;
    else
        error('Logic error!');
    end
elseif L<0
    if abs(L)<e
        L=0;
    else
        error('Logic error!');
    end
end
% ang=asin(L)/dis(u,[0;0])/dis(v,[0;0]);
ang=asin(L);
ang=ang/pi*180; % change to degrees
if u'*v>=0
    ;
else
    ang=180-ang; % This is due to asin() could only return a angle between 0 and 90 degrees
end