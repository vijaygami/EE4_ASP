clf; clc;
N=100;      % Number of data points
K=100;      % FFT size
FS=1000;    % Sample freq

sin=sine_gen(24, FS, N);    % Generate sinewave


subplot (1,3,1)
hold on
plot([1:N],sin,'.', 'color', [0 0.4470 0.7410])
plot([N+1: 2*N],sin, '.', 'color', [0 0.4470 0.7410])
plot(N, sin(N), '*', 'color', [0.8500 0.3250 0.0980])
plot(N+1, sin(1), '*', 'color', [0.8500 0.3250 0.0980])
hold off

title({'Discontinuity of 24Hz Sinewave', 
    '$\quad$ $f_s$=1000, N=100, K=100'}, 'FontWeight', 'normal')
xlabel('Sample Number')
ylabel('Sin(0.048$\pi$x)')
ylim([-1.2, 1.2])
box off; grid on


%_Increace K________________________________

N=100;      % Number of data points
K=1000000;  % FFT size
FS=1000;    % Sample freq
sin=sine_gen(24, FS, N);    % Generate sinewave

subplot (1,3,2)
hold on;
plot(x_axis(100, FS/2), abs(fftshift(fft(sin, N))), 'x');
plot(x_axis(K, FS/2), abs(fftshift(fft(sin, K))));
hold off;
title({'DFT, DTFT of 24Hz Sinewave', 
    '$\quad$ $f_s$=1000, N=100, K=100'}, 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Magnitude')
ylim([0,65]); xlim([-100, 100])
legend('K point DFT', 'DTFT', 'location', 'best')
box off; grid on

%_Coherent Sampling________________________________

N=125;      % Number of data points
K=125;      % FFT size
FS=1000;     % Sample freq
sin_coherent=sine_gen(24, FS, N);

subplot (1,3,3); hold on
plot(x_axis(K, FS/2), abs(fftshift(fft(sin_coherent, K))), 'x');
plot(x_axis(K*100, FS/2), abs(fftshift(fft(sin_coherent, K*100))));
title({'DFT, DTFT of 24Hz Sinewave (Coherent Sampling)',
    '$~~~~~~~~~~~~~~~~$ $f_s$=1000, N=125, K=125'}, 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Magnitude')
ylim([0,80]); xlim([-100, 100]);
legend('K point DFT', 'DTFT', 'location', 'best')

box off; grid on



