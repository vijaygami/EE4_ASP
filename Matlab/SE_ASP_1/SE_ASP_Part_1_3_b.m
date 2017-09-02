clc;clear;

N=256;  % Length of signal
L=16384; % length of FFT
variance=0;  % Variance of AWGN
a1=1;
a2=1;
phi1=0;
phi2=0;
f0=0.2;
alpha_set=[0.6, 0.65, 0.7, 0.75, 0.8, 0.9];

% periodogram for different alpha values
i=0;
figure(1); clf
for alpha = alpha_set
    
    i=i+1;
    subplot(1, size(alpha_set, 2), i);
    
    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    f=fft(x, L);
    pxx=(f.*conj(f))./L;
    pxx=fftshift(pxx);
    
    plot(x_axis(L, 0.5), 10*log10(pxx));
    xlabel({'Normalised Frequency',
    '~~~~~[x2$\pi$ rad/samp]'});
    ylabel('PSD [dB/rad/samp]');
    title(['Periodogram, $\alpha$ = ',num2str(alpha)], 'FontWeight', 'normal');
    xlim([0.187,0.215 ]);
    ylim([-17.5,1]);
    box off; grid on;
end


% Repeat with Hamming Window
figure(2); clf; i=0;
for alpha = alpha_set
    
    i=i+1;
    subplot(1, size(alpha_set, 2), i);

    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    xw=x.*hamming(N)';
    f=fft(xw, L);
    pxx=(f.*conj(f))./L;
    pxx=fftshift(pxx);
    
    plot(x_axis(L, 0.5), 10*log10(pxx));
    xlabel({'Normalised Frequency',
    '~~~~~[x2$\pi$ rad/samp]'});
    ylabel('PSD[dB/rad/samp]');
    title(['Periodogram, $\alpha$ = ',num2str(alpha)], 'FontWeight', 'normal');
    xlim([0.187,0.215 ]);
    ylim([-53,-4]);
    box off; grid on

end