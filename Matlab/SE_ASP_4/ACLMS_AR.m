function [ y_est, h, g, e ] = ACLMS_AR(y, u, p )
% Augmented Complex LMS that approximates the Wiener solution in a recursive fashion.
%
%Inputs:
%
% y = measured output of system to model. length N
% u = adaption gain
% p = adaptive filter order
%
%Output:
%
% y_est = LMS estimate of length N
% h = (p) * N matrix containing evolution of adaptive weights over time
% g = (p) * N matrix containing evolution of adaptive weights over time
% e = error vector e[n] = y[n] - y_est[n] 

N=length(y);        % length of input vector.
h=zeros(p,N);       % starting point for estimate
g=zeros(p,N);       % starting point for estimate
e=zeros(1,N);       % prealocate memory for speed
y_est=zeros(1,N);

for n=1:N-1,
    
    xp=xpast(y,n-1,p-1);                           % xp = (p) past values of x starting from x[n-1]. (column vector)
    y_est(1,n)= h(:,n)'*xp + g(:,n)'*conj(xp);     % estimated output

    e(1,n)=y(n)-y_est(1,n);                        % error in estimation of output
    h(:,n+1) = h(:,n) + u*conj(e(n))*xp;           % update filter weights
    g(:,n+1) = g(:,n) + u*conj(e(n))*conj(xp);     % update filter weights

end

h=h(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.
g=g(:, 1:N);

end

