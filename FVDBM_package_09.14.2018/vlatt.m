function [V,Mew]=vlatt(d,q,Tau)
% V=vlatt(d,q) defines the velocity decomposition of certain lattice
% d is the dimension of the domain
% q is the total number of velocity components during the decomposition
% Tau is the relaxation time
% V is a matrix contains x and y components of each velocity component
% Mew is the kinematic viscosity

if d==1
    ;
elseif d==2
    if q==5
        ;
    elseif q==7
        ;
    elseif q==9
        ;
    elseif q==13
        ;
    elseif q==37
        ;
    else
        error('Check lattice type!');
    end
    
    if q==5
        V=zeros(d,q);
        for l=1:q
            if l<=1
                V(:,l)=[0;0]; % Zero velocity of first component
            else
                V(:,l)=double(int16([cos(pi/2*(l-2));sin(pi/2*(l-2))]));
            end
        end
        Mew=Inf; % Still unknown
    end
    
    
    %%%% D2Q7
    if q==7
        V=zeros(d,q);
        for l=1:q
            if l<=1
                V(:,l)=[0;0]; % Zero velocity of first component
            else
                V(:,l)=[cos(pi/3*(l-2));sin(pi/3*(l-2))];
            end
        end
        Mew=Tau/4; 
    end
    
    
    %%%% D2Q9
    if q==9
        V=zeros(d,q);
        for l=1:q
            if l<=1
                V(:,l)=[0;0]; % Zero velocity of first component
            elseif l<=5
                V(:,l)=double(int16([cos(pi/2*(l-2));sin(pi/2*(l-2))]));
            else
                V(:,l)=double(int16(sqrt(2)*[cos(pi/2*(l-6)+pi/4);sin(pi/2*(l-6)+pi/4)]));
            end
        end
        Mew=Tau/3;
    end
    
    
    %%%% D2Q13
    if q==13
        V=zeros(d,q);
        for l=1:q
            if l<=1
                V(:,l)=[0;0]; % Zero velocity of first component
            elseif l<=5
                V(:,l)=double(int16([cos(pi/2*(l-2));sin(pi/2*(l-2))]));
            elseif l<=9
                V(:,l)=double(int16(sqrt(2)*[cos(pi/2*(l-6)+pi/4);sin(pi/2*(l-6)+pi/4)]));
            else
                V(:,l)=double(int16(2*[cos(pi/2*(l-10));sin(pi/2*(l-10))]));
            end
        end
        Mew=Tau/2;
    end
    
    
    %%%% D2Q37
    if q==37
        c_l=1.19697977039307435897239;
        V=zeros(d,q);
        for l=1:q
            if l<=1
                V(:,l)=[0;0]; % Zero velocity of first component
            elseif l<=5
                V(:,l)=double(int16([cos(pi/2*(l-2));sin(pi/2*(l-2))]))*c_l;
            elseif l<=9
                V(:,l)=double(int16(sqrt(2)*[cos(pi/2*(l-6)+pi/4);sin(pi/2*(l-6)+pi/4)]))*c_l;
            elseif l<=13
                V(:,l)=double(int16(2*[cos(pi/2*(l-10));sin(pi/2*(l-10))]))*c_l;
            elseif l<=17
                V(:,l)=double(int16(2*sqrt(2)*[cos(pi/2*(l-14)+pi/4);sin(pi/2*(l-6)+pi/4)]))*c_l;
            elseif l<=21
                V(:,l)=double(int16(3*[cos(pi/2*(l-18));sin(pi/2*(l-18))]))*c_l;
            else
                V(:,22)=[2;1]*c_l;
                V(:,23)=[1;2]*c_l;
                V(:,24)=[-1;2]*c_l;
                V(:,25)=[-2;1]*c_l;
                V(:,26)=[-2;-1]*c_l;
                V(:,27)=[-1;-2]*c_l;
                V(:,28)=[1;-2]*c_l;
                V(:,29)=[2;-1]*c_l;
                
                V(:,30)=[3;1]*c_l;
                V(:,31)=[1;3]*c_l;
                V(:,32)=[-1;3]*c_l;
                V(:,33)=[-3;1]*c_l;
                V(:,34)=[-3;-1]*c_l;
                V(:,35)=[-1;-3]*c_l;
                V(:,36)=[1;-3]*c_l;
                V(:,37)=[3;-1]*c_l;
            end
        end
        Mew=Tau;
    end
    
    
elseif d==3
    ;
else
    error('wrong dimension');
end