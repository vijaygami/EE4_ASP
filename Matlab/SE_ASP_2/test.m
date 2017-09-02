clc;clear;rng(1);
N=1000;
% 
% M=100
% 
% x=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
% [X,R] = corrmtx(x,M,'mod');
% 
% 
% x2=xcorr(x, M, 'unbiased');
% X2=toeplitz(x2(M+1:end));
% 
% 
% RR=X'*X;
% 
% e=R-RR;
% mean(mean(e))


w=bartlett(1025);

plot(real(fft(w)))
xlim([0,100])
