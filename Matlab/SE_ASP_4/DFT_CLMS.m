function [ y_est, w , e] = DFT_CLMS(y, N, gamma, mu)
% Adaptive DFT using leaky CLMS.
%
% Inputs:
%
% y = input vector with lenght of M
% N = DFT length (number of weigh vectors)
% gamma = leakage factor
% mu = adaptation gain, use 1 to converge to DFT solution in N steps
%
% Output:
%
% y_est = estimated value of y based on current weights 
% w = N * M matrix containing evolution of adaptive DFT weights over time
% e = error vector e[n] = y[n] - y_est[n] 
y=y(:);
M = length(y);


w = zeros(N,M);       % initialise weights and preallocate memory for speed
y_est = zeros(M, 1);  
e = zeros(M, 1);

for n = 0:M-1
    
    x(:,1) = (1/N)*(exp(1i*2*pi*n*[0:N-1]/N));
        
    y_est(n+1) = w(:, n+1)'*x;
    e(n+1) = y(n+1) - y_est(n+1);
    w(:, n+2) = (1-mu*gamma)*w(:,n+1)+mu*conj(e(n+1))*x;    % leaky lms weight update

end

w=w([1:N], [1:M]);

end


