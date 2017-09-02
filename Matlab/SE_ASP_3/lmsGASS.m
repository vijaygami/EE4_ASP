function [ y_est, w, e, u ] = lmsGASS( x, y, q, u0,  rho, alpha, sw)
% Adaptive filter that approximates the Wiener solution in a recursive
% fashion using GASS algorithms
%
% Inputs:
%
% x = input vector with lenght of N
% y = realistic (noisy) measured output of system to model. length N
% q = adaptive filter order
% u0 = initial adaption gain
% rho = step size for GASS algorithms
% alpha = Ang % Farhang GASS algo parameter. between 0 and 1.
% sw = switch for choosing the three GASS algos:
%
%       'A&F' for Ang & Farhang
%       'Ben' for Benveniste
%       'M&E' for Matthews and Xie
%
% Outputs:
%
% y_est = LMS estimate of length N
% e = error vector [n] = z[n] - y[n] 
% w = p*N matrix containing evolution of adaptive weights over time
% u = adaptation gains over time


N=length(x);              % length of input vector.

w=zeros(q+1,N);         % starting point for estimate
u=ones(1,N)*u0;         % initial adaptation gain 
phi=zeros(q+1,N);       % initial value and prealocate for speed
y_est=zeros(1,N);
e=zeros(1,N);

for n=q+1:N,
    
    xp=xpast(x,n,q);        % xp = (q) past values of x starting from x[n] = current value. (column vector)
    xp2=xpast(x,n-1,q);     % xp2 = (q) past values of x starting from x[n-1] = current value. (column vector)
    y_est(1,n)= xp'*w(:,n); % estimated output
    e(1,n)=y(n)-y_est(n);   % error in estimation of output

    switch sw
        case 'Ben'
            phi(:,n)=(eye(q+1)-u(n-1)*(xp2*xp2'))*phi(:,n-1) + e(n-1)*xp2;
        case 'A&F'
            phi(:,n)=alpha*phi(:,n-1) + e(n-1)*xp2;
        case 'M&X'
            phi(:,n)=e(n)*xp2;
        otherwise
            error('No valid GASS algo specified');
    end
    
    w(:,n+1) = w(:,n) + u(n)*e(n)*xp;  % update filter weights
    u(n+1)=u(n)+rho*e(n)*xp'*phi(:,n); % GASS algo     
    
end

% remove last update just to keep all vectors of length N.
w=w(:, 1:N);   
u=u(:, 1:N); 

end

