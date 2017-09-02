clc;clear;clf
%% FT
f=20;                       % frequency in Hz
t=[-0.01:0.0001:0.1];       % time vector
gt=sin(t*2*pi*f);           % sinewave


figure(1); clf; subplot(1,4,1)
plot(t, gt)
xlabel('Time [Seconds]')
ylabel('g(t)')
title('Infinite Sinewave', 'fontWeight', 'normal')
xlim([-0.01, 0.1]); ylim([-1.1, 1.1]); box off; grid on;

subplot(1,4,2); hold on
stem([-20, 20], [0.5, 0.5], '^', 'filled')
xlim([-40, 40]); ylim([0, 1])
hold off

xlabel('Frequency [Hz]')
ylabel('$|G(f)|$')
title('FT of Infinite Sinewave', 'fontWeight', 'normal')
box off; grid on;


%% DTFT
f=20;                       % frequency in Hz
fs=500;                     % sample frequency
t=[-0.01:1/fs:0.1];         % time vector
gt=sin(t*2*pi*f);           % sinewave
windowed=[zeros(1,5), gt(6:40), zeros(1,16)];

subplot(1,4,3); hold on
stem(t, windowed)
plot([0,0], [0,1], 'color', 'r')
plot([0.068,0.068], [0,1], 'color', 'r')
plot([0,0.068], [1,1], 'color', 'r')
plot([0.068, 0.1], [0,0], 'color', 'r')
plot([-0.1, 0], [0,0], 'color', 'r')
hold off

legend('Windowed Sinewave', 'Window ($T_w$ Samples Long)', 'location', 'best')
xlabel('Time [Seconds]')
ylabel('g(t)')
title('Discrete Windowed Sinewave', 'fontWeight', 'normal')
xlim([-0.01, 0.1]); ylim([-1.1, 1.7]); box off; grid on;

L=4096;
subplot(1,4,4);
plot(x_axis(L, fs/2), fftshift(abs(fft(gt(6:40), L))))
xlim([-fs/2, fs/2]); ylim([0, 20])
box off; grid on

xlabel('Frequency [Hz]')
ylabel('$|G(f)|$')
title('DTFT of Windowed Sinewave', 'fontWeight', 'normal')

