function [ y, w, e ] = lms( x, z, u, q )
% modified from my (Vijay Gami) ASP courcework year 3
% Adaptive filter that approximates the Wiener solution in a recursive
% fashion.
%
%Inputs:
%
% x = input vector with lenght of N
% z = measured output of system to model (may have additive noise). length N
% u = adaption gain
% q = adaptive filter order
%
%Output:
%
% y = LMS estimate of length N
% e = error vector [n] = z[n] - y[n] 
% w = p*N matrix containing evolution of adaptive weights over time


N=length(x);        % length of input vector.
w=zeros(q+1,1);     % starting point for estimate
e=zeros(1,N);       % prealocate memory for speed
y=zeros(1,N);

for n=1:N,
   
    xp=xpast(x,n,q);    % xp = (q) past values of x starting from x[n] = current value. (column vector)
    y(1,n)= xp'*w(:,n); % estimated output
    e(1,n)=z(n)-y(1,n); % error in estimation of output
    w(:,n+1) = w(:,n) + u*e(n)*xp; % update filter weights

end

w=w(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.

end

