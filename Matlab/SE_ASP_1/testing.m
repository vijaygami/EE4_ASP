clc;clear
rng(1)
N=2000;
A=2000;
x=randn(A, N);
var=5;
x=x*sqrt(5);



for i=1:A

[p(i,:), w] = periodogram(x(i,:), rectwin(N), N, 'twosided');

pp(i,:)=abs(fftshift(pgm(x(i,:), N)));

end

mean(mean(p))
mean(mean(pp))


mean(std(p).^2)
mean(std(pp).^2)
