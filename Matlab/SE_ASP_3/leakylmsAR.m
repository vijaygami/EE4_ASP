function [ y, e, a ] = leakylmsAR( x, u, p, gamma )
% Code modified from my (Vijay Gami) code from ASP year 3 CW.
% Adaptive filter that adaptivley identifies AR system in a recursive
% fashion using the leaky LMS algorithm
%
% Inputs:
%
% x = input vector with lenght of N
% u = adaption gain
% p = adaptive filter order .(number of co-eff is order+1)
% gamma = leakage coefficient. (between 0 and 1)
%
% Outputs:
%
% y = predictor output
% e = error vector [n] = x[n] - y[n] 
% a = (p) * N matrix containing evolution of adaptive weights over time,
%  it does not include a0 which is 1.

N=length(x);        % length of input vector.
a=zeros(p,1);       % starting point for estimate
e=zeros(1,N);       % prealocate memory for speed
y=zeros(1,N);

for n = 1:N,
    
    xp = xpast(x,n-1,p-1);          % xp = (p) past outputs of x, NOT including  x[n], hence i-1, p-1 inputted. (column vector)
    y(1,n) = xp'*a(:,n);            % estimated output
    e(1,n) = x(n)-y(1,n);           % error in estimation of output
    a(:,n+1) = (1-u*gamma)*a(:,n) + u*e(1,n)*xp;  % update all filter weights excluding the first weight. 

end

a=a(:, 1:N);   % remove last coefficient update just to keep all vectors of length N.

end

