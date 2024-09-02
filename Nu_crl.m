function Nu=Nu_crl(Ra,flag)

%% This function calculates the Nu based on Ra and the type of correlation

if flag==0 %% Experiment, Hollands et. al (1976), Air (Pr=0.71) Ra<=10e+8
    if Ra>1e+8
        error('Exceed the limit! Choose a Ra<1e+8');
    end
    a=1-1708/Ra;
    b=Ra^(1/3)/18-1;
    if a<=0
        a=0;
    end
    if b<=0
        b=0;
    end
    Nu=1+1.44*a+b;
elseif flag==1 %% from the textbook, 0.5<Pr<2, 
    if Ra>3.2e+5 || Ra<1700
        error('Exceed the bounds! Choose a Ra<3.2e+5 and Ra>1700');
    end
    if Ra<=7000
        Nu=0.059*Ra^0.4;
    else
        Nu=0.212*Ra^(1/4);
    end
else
    error('Wrong flag for Nu correlations!');
end