clc;clear; clf
M=10;
L=256;

x1=(1/M)*[[1:1:M-1],[M:-1:1],zeros(1,(L-(2*M-1)))];  % Zero padded vector of length L
x1f=fftshift(fft(x1));

M=128;
x2=(1/M)*[[1:1:M-1],[M:-1:1],zeros(1,(L-(2*M-1)))];  % Zero padded vector of length L
x2f=fftshift(fft(x2));

subplot(2,2,1)
plot(x1)
title('Incorrectly Shifted ACF, M=10, L=256','FontWeight', 'normal')
xlabel('Sample Number')
ylabel('z')
box off; grid on


subplot(2,2,3)
plot(x_axis(L, 0.5) , real(x1f))
title('Real(DFT(z)), M=10, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [cycles/samp]')
ylabel('Real(DFT(z)))')
box off; grid on
xlim([-0.5, 0.5])
ax = gca;
ax.XTick = [-0.5:0.125:0.5];
subplot(2,2,2)
plot(x2)
title('Incorrectly Shifted ACF, M=128, L=256','FontWeight', 'normal')
xlabel('Sample Number')
ylabel('z')
box off; grid on

subplot(2,2,4)
plot(x_axis(L, 0.5) , real(x2f))
title('Real(DFT(z)), M=128, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [cycles/samp]')
ylabel('Real(DFT(z)))')
box off; grid on
ylim([-60,150])
xlim([-0.5, 0.5])
ax = gca;
ax.XTick = [-0.5:0.125:0.5];