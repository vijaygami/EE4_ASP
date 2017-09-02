function [ y_est, h, e ] = CLMS(x, y, u, q )
% Complex LMS that approximates the complex Wiener solution in a recursive fashion.
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
% y = LMS estimate of length N
% h = (q+1) * N matrix containing evolution of (conjugate) adaptive weights over time
% e = error vector e[n] = y[n] - y_est[n] 


N=length(x);        % length of input vector.
h=zeros(q+1,N);     % starting point for estimate & prealocate memory for speed
e=zeros(1,N);       
y_est=zeros(1,N);

for n=1:N,
    
    xp = xpast(x,n,q);                    % xp = (q) past values of x starting from x[n] = current value. (column vector)
    y_est(1,n) = h(:,n)'*xp;              % estimated output

    e(1,n)=y(n)-y_est(1,n);               % error in estimation of output
    h(:,n+1) = h(:,n) + u*conj(e(n))*xp;  % update filter weights

end

h=h(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.

end

