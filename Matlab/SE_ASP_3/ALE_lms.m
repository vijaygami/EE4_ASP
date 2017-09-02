function [ x_est, e, w ] = ALE_lms( s, delta, mu, M )
% ALE using LMS for the linear prediction
%
%Inputs:
%
% s = noisy signal
% delta = delay for decorelation of noise
% mu = adaption gain of LMS algo
% M = adaptive filter order
%
%Output:
%
% x_est = estimated de-noised signal
% e = error  e[n] = s[n] - x_est[n] 
% w = (p+1) * N matrix containing evolution of adaptive weights over time


N=length(s);            % Length of input vector.
w=zeros(M+1,1);         % Starting point for estimate
x_est=zeros(1,N);       % Preaclocate memory for speed
e=zeros(1,N);

delayed = [zeros(1,delta), s];

for n=1:N,
    
    up=xpast(delayed,n,M);          % up = (M) past values of x starting from x[n] = current value. (column vector)
    x_est(:,n)= up'*w(:,n);         % estimated output
    e(1,n)=s(n)-x_est(1,n);         % error in estimation of output
    w(:,n+1) = w(:,n) + mu*e(n)*up; % update filter weights
end

w=w(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.
end

