function y=sym_para(x1,x2,y_max,x)
% function y=sym_para(x1,x2,y_max,x)
% Generate the y value at any location x
% from the equation y=a*x^2+b*x+c that satisfies
% y=0 at x=x1 and x2 and y=y_max at x=(x1+x2)/2
% So the profile of y is symmetric between x1 and x2 with y_max the maximum
% value

% x1 is the lower bound of independent varible
% x2 is the upper bound of independent varible
% y_max is the maximum value of y
% x is the any location between x1 and x2
% y is the generated value at x

if length(x1)>1 || length(x2)>1 || length(y_max)>1 || length(x)>1
    error('Check the dimension of the inputs!');
end

if x1>=x2
    error('The lower bound cannot be bigger than the upper bound!');
end

if x<x1 || x>x2
    error('The location is out of range!');
end

a=-4*y_max/(x1-x2)^2;
b=4*y_max*(x1+x2)/(x1-x2)^2;
c=-4*y_max*x1*x2/(x1-x2)^2;

y=a*x^2+b*x+c;

e=(x1+x2)/2*10;
if single(y+e)==single(e)
    y=0;
end