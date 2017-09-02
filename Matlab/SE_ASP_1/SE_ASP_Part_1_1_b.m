clc;clear;clf
N=100;      % Number of data points
K=100;      % FFT size
FS=1000;    % Sample freq

sin=sine_gen(20, FS, N);    % Generate sinewave

subplot(1,2,1)
stem(x_axis(K, FS/2), (1/K)*abs(fftshift(fft(sin, K))));
title('(Scaled) DFT of 20Hz Sinewave, $f_s$=1000, N=100, K=100', 'FontWeight', 'normal')
xlabel('Frequency Bin [Hz]')
ylabel('Magnitude')
ylim([0,.6])
grid on; box off

K=1000;      % FFT size
subplot(1,2,2);
stem(x_axis(K, FS/2), abs(fftshift(fft(sin, K))));
title('DFT of 20Hz Sinewave, $f_s$=1000, N=100, K=1000', 'FontWeight', 'normal')
xlabel('Frequency Bin [Hz]')
ylabel('Magnitude')
xlim([-80, 80])
ylim([0,60])
grid on; box off