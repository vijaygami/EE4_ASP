clc;clear;
rng(2);
L=128;              % FFT length (zero pad)

figure(1); clf;
%% __0.3 and 0.32 hZ sinusoids__________
N=30;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;

subplot(2,3,1)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.32']}, 'FontWeight', 'normal')

N=35;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;

subplot(2,3,2)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.32']}, 'FontWeight', 'normal')

N=45;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;

subplot(2,3,3)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.32']}, 'FontWeight', 'normal')


%% __N=30__________
rng(2)
N=25;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;

subplot(2,3,4)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.32']}, 'FontWeight', 'normal')

N=25;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.33*n)+noise;

subplot(2,3,5)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.33']}, 'FontWeight', 'normal')

N=25;
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.34*n)+noise;

subplot(2,3,6)
plot(x_axis(L, 0.5)+0.5, 10*log10(fftshift(pgm(x, L))))
title({'~~~~~~~~~~~~~~~~Periodogram:', 
    ['N = ', num2str(N), ', L = ', num2str(L), ', $f_0$ = 0.3, $f_1$ = 0.34']}, 'FontWeight', 'normal')


% setting graph axis and settings in a compact manor
for i=1:6
    subplot(2,3,i)
    box off
    xlabel('Frequency [Hz]')
    ylabel('PSD [dB/Hz]')
    ylim([-31, 18])
    grid on

end