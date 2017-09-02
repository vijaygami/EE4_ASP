function [ y_est, h, e ] = CLMS_AR(y, u, p)
% Complex LMS
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
% h = p * N matrix containing evolution of adaptive weights over time
% e = error vector e[n] = y[n] - y_est[n] 

N=length(y);        % length of input vector.
h=zeros(p,N);       % starting point for estimate & prealocate memory for speed
e=zeros(1,N);       
y_est=zeros(1,N);

for n=1:N,
    
    xp = xpast(y,n-1,p-1);                % xp = (p) past values of x starting from x[n-1] (column vector)
    y_est(1,n) = xp.'*conj(h(:,n));       % estimated output

    e(1,n)=y(n)-y_est(1,n);               % error in estimation of output
    h(:,n+1) = h(:,n) + u*conj(e(n))*xp;  % update filter weights

end

h=h(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.

end

