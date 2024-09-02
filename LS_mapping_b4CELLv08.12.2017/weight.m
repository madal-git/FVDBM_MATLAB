function [wh,wt]=weight(d,qh,qt)
% [wh,wt]=weight(d,qh,qt) return the weighting vector
% for calculating equilibrium pdf.
% d is the dimension of the domain
% qh is the total number of hydrodynamic pdf components during 
% the decomposition
% qt is the total number of thermal pdf components during 
% the decomposition
% wh is the corresponding hydrodynamic weighting vector
% wt is the corresponding thermal weighting vector

%%%%%% Hydrodynamic
if d==2
    if qh==7 % D2Q7
        wh(1)=1/2;wh(2)=1/12;wh(3)=1/12;wh(4)=1/12;wh(5)=1/12;wh(6)=1/12;wh(7)=1/12;
    elseif qh==9 % D2Q9
        wh(1)=4/9;wh(2)=1/9;wh(3)=1/9;wh(4)=1/9;wh(5)=1/9;wh(6)=1/36;wh(7)=1/36;wh(8)=1/36;wh(9)=1/36;
    elseif qh==13 % D2Q13
        wh(1)=3/8;wh(2)=1/12;wh(3)=1/12;wh(4)=1/12;wh(5)=1/12;wh(6)=1/16;wh(7)=1/16;wh(8)=1/16;wh(9)=1/16;wh(10)=1/96;wh(11)=1/96;wh(12)=1/96;wh(13)=1/96;
    elseif qh==37
        wh(1)=0.2331506691323525022865;
        
        wh(2)=0.10730609154221900241246;
        wh(3)=0.10730609154221900241246;
        wh(4)=0.10730609154221900241246;
        wh(5)=0.10730609154221900241246;
        
        wh(6)=0.05766785988879488203006;
        wh(7)=0.05766785988879488203006;
        wh(8)=0.05766785988879488203006;
        wh(9)=0.05766785988879488203006;
        
        wh(10)=0.01420821615845075026469;
        wh(11)=0.01420821615845075026469;
        wh(12)=0.01420821615845075026469;
        wh(13)=0.01420821615845075026469;
        
        wh(14)=0.0010119375926735754754;
        wh(15)=0.0010119375926735754754;
        wh(16)=0.0010119375926735754754;
        wh(17)=0.0010119375926735754754;
        
        wh(18)=0.00024530102775771734547;
        wh(19)=0.00024530102775771734547;
        wh(20)=0.00024530102775771734547;
        wh(21)=0.00024530102775771734547;
        
        wh(22)=0.00535304900051377523273;
        wh(23)=0.00535304900051377523273;
        wh(24)=0.00535304900051377523273;
        wh(25)=0.00535304900051377523273;
        wh(26)=0.00535304900051377523273;
        wh(27)=0.00535304900051377523273;
        wh(28)=0.00535304900051377523273;
        wh(29)=0.00535304900051377523273;
        
        wh(30)=0.0002834142529941982174;
        wh(31)=0.0002834142529941982174;
        wh(32)=0.0002834142529941982174;
        wh(33)=0.0002834142529941982174;
        wh(34)=0.0002834142529941982174;
        wh(35)=0.0002834142529941982174;
        wh(36)=0.0002834142529941982174;
        wh(37)=0.0002834142529941982174;
    else
        error('Temporarily not available!');
    end
elseif d==3
    error('Temporarily not available!');
else
    error('Dimension is wrong!');
end
%%%%%% Thermal for Passive-Scalar
if d==2
    if qt==7 % D2Q7
        error('Temporarily not available!');
    elseif qt==9 % D2Q9
        wt(1)=-2/3;wt(2)=1/9;wt(3)=1/9;wt(4)=1/9;wt(5)=1/9;wt(6)=1/36;wt(7)=1/36;wt(8)=1/36;wt(9)=1/36;
    elseif qt==13 % D2Q13
        error('Temporarily not available!');
    else
        error('Wrong lattice type for thermal model!');
    end
elseif d==3
    error('Temporarily not available!');
else
    error('Dimension is wrong!');
end