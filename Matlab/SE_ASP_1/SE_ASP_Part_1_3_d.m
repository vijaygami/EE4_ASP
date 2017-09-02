clc;clear;

N=256;       % Length of signal
L=16384;     % length of FFT
variance=0;  % Variance of AWGN
a1=1;
a2_set=[1, 0.1, 0.01, 0.001];   % Set of values of a2 to plot for
phi1=0;
phi2=0;
f0=0.2;
alpha=4;


% periodogram for different a2 values
figure(1); clf; i=0;
for a2 = a2_set
    
    i=i+1;
    subplot(1, size(a2_set, 2), i);
    
    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    f=fft(x, L);
    pxx=(f.*conj(f))./L;
    pxx=fftshift(pxx);
     
    plot(x_axis(L, 0.5), 10*log10(pxx));
    xlabel({'Normalised Frequency',
    '~~~~~[x2$\pi$ rad/samp]'});
    ylabel('PSD [dB/rad/samp]');
    title(['Periodogram, a2 = ',num2str(a2), ', $\alpha$ = ', num2str(alpha)], 'FontWeight', 'normal');
    xlim([0.192,0.25 ]);
    ylim([-33, 2]);
    box off
    grid on
    
end


% Repeat with alpha = 12
figure(2); clf; i=0;
alpha=12;
for a2 = a2_set
    
    i=i+1;
    subplot(1, size(a2_set, 2), i);
    
    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    f=fft(x, L);
    pxx=(f.*conj(f))./L;
    pxx=fftshift(pxx);
    
    plot(x_axis(L, 0.5), 10*log10(pxx));
    xlabel({'Normalised Frequency',
    '~~~~~[x2$\pi$ rad/samp]'});
    ylabel('PSD [dB/rad/samp]');
    title(['Periodogram, a2 = ',num2str(a2), ', $\alpha$ = ', num2str(alpha)], 'FontWeight', 'normal');
    xlim([0.192,0.27 ]);
    ylim([-35, 2]);
    box off
    grid on
    
end

