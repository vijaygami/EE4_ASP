function [ y_est, h, g, e ] = ACLMS(x, y, u, q )
% Augmented Complex LMS.
%
%Inputs:
%
% x = input vector with lenght of N
% y = measured output of system to model (may have additive noise). length N
% u = adaption gain
% q = adaptive filter order
%
%Output:
%
% y_est = LMS estimate of length N
% h = (q+1) * N matrix containing evolution of adaptive weights over time
% g = (q+1) * N matrix containing evolution of adaptive weights over time
% e = error vector e[n] = y[n] - y_est[n] 

N=length(x);        % length of input vector.
h=zeros(q+1,N);     % starting point for estimate
g=zeros(q+1,N);     % starting point for estimate
e=zeros(1,N);       % prealocate memory for speed
y_est=zeros(1,N);

for n=1:N-1,
    
    xp=xpast(x,n,q);                               % xp = (q) past values of x starting from x[n] = current value. (column vector)
    y_est(1,n)= h(:,n)'*xp + g(:,n)'*conj(xp);     % estimated output

    e(1,n)=y(n)-y_est(1,n);                        % error in estimation of output
    h(:,n+1) = h(:,n) + u*conj(e(n))*xp;           % update filter weights
    g(:,n+1) = g(:,n) + u*conj(e(n))*conj(xp);     % update filter weights

end

h=h(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.
g=g(:, 1:N);
end

