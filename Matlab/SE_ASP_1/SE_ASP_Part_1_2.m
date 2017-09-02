clc;clear;clf;
N=100;       % Number of data points
fs=1000;     % Sample freq

sig=sine_gen(100, fs, N);    % Generate sinewave
%sig=[zeros(1, N/2),ones(1,N/2)];
%sig=randn(1,N);

figure(1)
% _PSD estimate with FFT(ACF) method
subplot(1,3,1)
xcor=xcorr(sig, 'unbiased');
xcor=ifftshift(xcor);
pxx_1=abs(fftshift(fft(xcor)));
plot(x_axis(length(xcor), fs/2), 10*log10(pxx_1))
title({'~~~~~~DFT of Unbiased ACF Estimate',
    '100Hz Sinewave, $f_s$=1000, N=100, K=100'}, 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
box off; grid on
ylim([-40, 20])

% _PSD estimate with FFT (biased ACF) method
subplot(1,3,2)
xcor=xcorr(sig, 'biased');
xcor=ifftshift(xcor);
pxx_2=abs(fftshift(fft(xcor)));
plot(x_axis(length(xcor), fs/2), 10*log10(pxx_2))
title({'~~~~~~~DFT of Biased ACF Estimate',
    '100Hz Sinewave, $f_s$=1000, N=100, K=100'}, 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
box off; grid on
ylim([-40, 20])



% _PSD estimate with direct FFT method 
subplot(1,3,3)
pxx_3=pgm(sig, 2*N-1);
plot(x_axis(length(pxx_3), fs/2), 10*log10(pxx_3))
title({'Periodogram of 100Hz Sinewave', 
    '~~~~~$f_s$=1000, N=100, K=100'}, 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
box off; grid on
ylim([-40, 20])




e=(pxx_1-pxx_3);

figure(2);
plot(e)

