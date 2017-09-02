% clc;clear
% % N=1000;
% % rng(1)
% % b=[1.5+1i, 2.5-0.5i];       % array of MA coef for sythesis stage
% % 
% % load('low-wind.mat'); v_h = v_east + 1i*v_north;
% %     cwgn = wgn(1,N,0,'complex');
% %     y=zeros(1,N);
% % 
% %     for k=2:N
% %         y(k) = cwgn(k) + b(1)*cwgn(k-1) + b(2)*conj(cwgn(k-1));
% %     end 
% % 
% % % ccoef=mean(y.*y)./mean(y.*conj(y))
% % % ccoef=mean(cwgn.*cwgn)./mean(cwgn.*conj(cwgn))
% % % ccoef=mean(v_h.*v_h)./mean(v_h.*conj(v_h))
% % 
% % cv=cov(y, conj(y));
% % cv=cv(1,2);
% % ccoef= cv/(var(y).*var(conj(y)))
% 
% 
% % F = [ 1, 1, 1; 1, exp(j*2*pi/3), exp(j*2*(2)*pi/3) ; 1, exp(j*2*2*pi/3), exp(j*2*2*2*pi/3)];
% % 
%  x=[1, 2, 3].';    % input
% % 
% % w=((F'*F)^-1)*F'*x;
% % ww=fft(x)./3;
% % 
% % w-ww
% % 
% % mat=((F'*F)^-1)*F';
% % mat2=F^-1;
% % 
% % mat-mat2;
% 
% F = [ 1, 1, 1; 1, exp(j*2*pi/3), exp(j*2*(2)*pi/3) ; 1, exp(j*2*2*pi/3), exp(j*2*2*2*pi/3)];
% 
% f=3*F^-1*x
% f-fft(x)

signal=f_clms(4,[300:end]);
signal=detrend(signal);

plot(x_axis(10000, 500), abs(fftshift(fft(signal, 10000))));
