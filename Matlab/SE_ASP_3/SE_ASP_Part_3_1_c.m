clc;clear
rng(1)

N=5000;             % number of input samples
p=2;                % AR process order;
u=[0.01 0.05];      % adaptation gain
o2=0.25;            % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
b1=[1];             % AR model for sythesis stage
a1=[1, -0.1, -0.8];

%% mean square error

n_ave=100;    % Number of averages to perform
MSE=zeros(2,N);  % Will store MSE

for j=1:n_ave  
    % synthesis stage
    rng(j)   
    n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input, zero mean and variance = o2
    x=filter(b1,a1,n); x=x(501:end, 1);                     % output of filter. Discard first 500 samples due to transient response

    for i = 1:length(u)
        [ y, e(i,:), a ] = lmsAR( x, u(i), p);   
        e(i,:)=e(i,:).^2;
        MSE(i,:)=MSE(i,:)+e(i,:);
    end
    
end

display('Experimental Time averaged steady state MSE')
MSE=MSE./n_ave;                         % Mean of error squared
TA_MSE=mean(MSE(:, [4000:end]), 2)      % Time average of steady state values only

display('Experimental Excess MSE')
MLMS = TA_MSE - o2

display('Experimental Missadjustment')
MLMS=MLMS./o2


display('Theoretical Missadjustment (Approx)')
% Theoretical misadjustment
% rxx=xcorr(x,1, 'unbiased')   %only need ACF from lags 0 to lag 1.
% Rxx=toeplitz(rxx(2:3))       %ACF matrix for lags 0 to 1 toeplitz matrix

Rxx = [25/27, 25/54 ; 25/54, 25/27]; % beter to use actual ACF determined in part 3.1.a rather than an estimate
MLMS_T_A = [0.5*u.*trace(Rxx)]'

display('Theoretical Missadjustment (Exact)')
lamda=eig(Rxx);
for i = 1:length(u)
    MLMS_T_E(i,1) = sum(u(i)*lamda./(2-u(i)*lamda));
end

display(MLMS_T_E)

