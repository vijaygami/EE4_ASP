function [ y_est, w, e ] = GNGD( x, y, q, u, rho, beta)
% Adaptive filter that approximates the Wiener solution in a recursive
% fashion.
%
%Inputs:
%
% x = input vector with lenght of N
% y = realistic (noisy) measured output of system to model. length N
% q = adaptive filter order .(number of co-eff is order+1)
% u = adaption gain
% rho = step for the adaption of regularisation factor (epsilon)
% beta = extra parameter fror scaling adaptation gain. Set to 1 for most cases.

%Output:
%
% y_est = LMS estimate of length N
% e = error vector e[n] = y[n] - y_est[n] 
% w = (p+1) * N matrix containing evolution of adaptive weights over time

 
N=length(x);              % length of input vector.
w=zeros(q+1,N);           % starting point for estimate
epsilon=ones(1,N)*u*rho;  % initial value and prealocate for speed
y_est=zeros(1,N);
e=zeros(1,N);

for n=q+1:N,
    
    xp=xpast(x,n,q);        % xp = (q) past values of x starting from x[n] = current value. (column vector)
    xp2=xpast(x,n-1,q);     % xp2 = (q) past values of x starting from x[n-1] = current value. (column vector)
    y_est(1,n)= xp'*w(:,n); % estimated output
    e(1,n)=y(n)-y_est(n);   % error in estimation of output
    
    w(:,n+1) = w(:,n) + (beta/(epsilon(n)+(xp'*xp)))*e(n)*xp; %update filter weights
    epsilon(n+1) = epsilon(n) - rho*u*((e(n)*e(n-1)*xp'*xp2)/((epsilon(n-1)+(xp2'*xp2))^2)); % update regularisation factor

end
w=w(:, 1:N); % remove last update just to keep all vectors of length N.

end

