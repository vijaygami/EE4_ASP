M=10;
L=256;

x1=(1/M)*[[M:-1:1],zeros(1,(L-(2*M-1))),[1:1:M-1]];  % Zero padded vector of length L
x1f=fft(x1);

M=128;
x2=(1/M)*[[M:-1:1],zeros(1,(L-(2*M-1))),[1:1:M-1]];  % Zero padded vector of length L
x2f=fft(x2);

% plots comparing real part vs magnitude
figure(1)
subplot(1, 4, 1);
plot(x_axis(length(x1f), pi),abs(x1f));
title('Abs (FFT(x)), M=10, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [rad/samp]')
ylabel('Amplitude')

subplot(1, 4, 2);
plot(x_axis(length(x1f), pi),real(x1f));
title('Real(FFT(x)), M=10, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [rad/samp]')
ylabel('Amplitude')


subplot(1, 4, 3);
plot(x_axis(length(x1f), pi),abs(x2f));
title('Abs (FFT(x)), M=128, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [rad/samp]')
ylabel('Amplitude')

subplot(1, 4, 4);
plot(x_axis(length(x1f), pi),real(x2f));
title('Real(FFT(x)), M=128, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [rad/samp]')
ylabel('Amplitude')


ratio=(imag(x1f)./real(x1f));
ratio(isnan(ratio))=0;          % NAN occurs when we have 0/0
mean(ratio)
var (ratio)

ratio=(imag(x2f)./real(x2f));
ratio(isnan(ratio))=0;
mean(ratio)
var (ratio)


max(imag(x1f)./real(x1f))
max(imag(x2f)./real(x2f))


% Plots of imaginary parts only
figure(2)
plot(x_axis(256, 0.5), imag(fftshift(x2f)))
title('Imaginary part of FFT(x), M=128, L=256','FontWeight', 'normal')
xlabel('Normalised Frequency [cycles/samp]')
ylabel('imag(FFT(x))')
xlim([-0.5, 0.5])
ylim([-3.5e-15, 3.5e-15])

box off; grid on
