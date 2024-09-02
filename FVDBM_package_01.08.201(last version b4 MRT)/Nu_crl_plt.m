function Nu_crl_plt(Ra_low,Ra_high,flag)

%% This function plot the Nu vs. Ra based on the used correlation
N=100;
Ra=Ra_low:(Ra_high-Ra_low)/(N-1):Ra_high;
Nu=zeros(1,N);
for i=1:N
    Nu(i)=Nu_crl(Ra(i),flag);
end
semilogx(Ra,Nu)
grid on;
xlabel('Ra')
ylabel('Nu')