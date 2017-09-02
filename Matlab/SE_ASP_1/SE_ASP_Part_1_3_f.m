clc;clear;

N=256;       % Length of signal
L=16384;     % length of FFT
variance=0;  % Variance of AWGN
a1=1;
%a2_set=[1, 0.1, 0.01, 0.001];   % Set of values of a2 to plot for
a2_set=[0.1, 0.01, 0.001, 0.0001];   % Set of values of a2 to plot for

phi1=0;
phi2=0;
f0=0.2;
alpha=4;
sidelobe=120; % Attenuation (in dB) of sidelobes of chebychev window compared to mainlobe
M=150;        % Window length for blackman tukey


% Periodogram using chebycshev window and coresponding blackman tukey method
i=0;
figure(1); clf

for a2 = a2_set
 
    i=i+1;

    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    
    % Chebychev window periodogram
    [pxx_c,w] = periodogram(x,chebwin(N, sidelobe), L, 1, 'twosided');
    pxx_c=fftshift(pxx_c);
    
    % Blackman Tukey method with Chebyshev window
    x_cor = xcorr(x, 'biased');                          % Biased ACF
    %x_cor = xcorr(x, 'unbiased');                       % Unbiased ACF
    x_cor = x_cor(N-M:N+M);                              % Use 2M+1 values of the ACF
    x_cor = x_cor.*chebwin(length(x_cor), sidelobe)';    % Chebychev window
    x_cor = ifftshift(x_cor);
    k=length(x_cor); k=(k+1)/2;
    x_cor = [x_cor(1:k), zeros(1,L-2*k+1), x_cor(k+1:2*k-1)];  % Zero padd ACF while maintaining symmetry
    pxx_bt=fftshift(fft(x_cor));                               % PSD estimate
    
    % Plot both estimates on same graph
    subplot(2, size(a2_set, 2), i);
    hold on
    plot(x_axis(length(pxx_c), 0.5), 10*log10(pxx_c),'LineWidth', 1)
    plot(x_axis(length(pxx_bt), 0.5), 10*log10(pxx_bt), 'LineWidth', 1)
    hold off
    
    xlabel('Normalised Frequency [cycles/samp]');
    ylabel('PSD [dB/cycles/samp]');
    title(['Periodogram: a2=',num2str(a2), ', $\alpha$=', num2str(alpha), ', ', num2str(sidelobe), 'dB'], 'FontWeight', 'normal');
    xlim([0.192,0.25 ]);
    ylim([-150, 55]);
    legend ('Chebyshev','Blackman-Tukey', 'Orientation','vertical', 'Location', 'best')
    grid on

end


% Another value of alpha
alpha = 12;
for a2 = a2_set
 
    i=i+1;

    x = a1*sin(f0*pi*2*(0:N-1) + phi1) + a2*sin((f0 + alpha/N)*pi*2*(0:N-1) + phi2) + sqrt(variance)*randn(1, N);
    
    % Chebychev window periodogram
    [pxx_c,w] = periodogram(x,chebwin(N, sidelobe), L, 1, 'twosided');
    pxx_c=fftshift(pxx_c);
    
    % Blackman Tukey method with Chebyshev window
    x_cor = xcorr(x, 'biased');                          % Biased ACF
    %x_cor = xcorr(x, 'unbiased');                       % Unbiased ACF
    x_cor = x_cor(N-M:N+M);                              % Use 2M+1 values of the ACF
    x_cor = x_cor.*chebwin(length(x_cor), sidelobe)';    % Chebychev window
    x_cor = ifftshift(x_cor);
    k=length(x_cor); k=(k+1)/2;
    x_cor = [x_cor(1:k), zeros(1,L-2*k+1), x_cor(k+1:2*k-1)];  % Zero padd ACF while maintaining symmetry
    pxx_bt=fftshift(fft(x_cor));                               % PSD estimate
    
    % Plot both estimates on same graph
    subplot(2, size(a2_set, 2), i);
    hold on
    plot(x_axis(length(pxx_c), 0.5), 10*log10(pxx_c),'LineWidth', 1)
    plot(x_axis(length(pxx_bt), 0.5), 10*log10(pxx_bt), 'LineWidth', 1)
    hold off
    
    xlabel('Normalised Frequency [cycles/samp]');
    ylabel('PSD [dB/cycles/samp]');
    title(['Periodogram: a2=',num2str(a2), ', $\alpha$=', num2str(alpha), ', ', num2str(sidelobe), 'dB'], 'FontWeight', 'normal');
    xlim([0.192,0.27 ]);
    ylim([-150, 50]);
    legend ('Chebyshev','Blackman-Tukey', 'Orientation','vertial', 'Location', 'best')
    grid on

end


