function [ x_est, w ] = ANC_lms( s, epsilon, mu, M )
% Adaptive Noise Cancelation using LMS adaptive filter
%
% Inputs:
% s = noisy signal of  lenght of N
% epsilon = secondary noise input length N
% mu = adaption gain
% M = adaptive filter order
%
% Output:
%
% x_est = estimated de-noised signal
% w = (M+1) * N matrix containing evolution of adaptive weights over time


N=length(s);            % Length of input vector.
w=zeros(M+1,1);         % Starting point for estimate
x_est=zeros(1,N);       % Preaclocate memory for speed
n_est=zeros(1,N);

for n=1:N,
    
    ep=xpast(epsilon,n,M);                  % xp = (M) past values of x starting from x[n] = current value. (column vector)
    n_est(1,n)= ep'*w(:,n);                 % estimated noise
    x_est(1,n)=s(n)-n_est(1,n);             % error in estimation of output
    w(:,n+1) = w(:,n) + mu*x_est(n)*ep;     % update filter weights

end

w=w(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.
end

