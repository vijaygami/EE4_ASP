clc;clear;clf;

N=128;        % Length of signal
variance=1;   % Variance of AWGN
a1=1.5;       % Amplitude of sinusoid
a2=0.75;
f0=0.15;      % (Normalised [Hz/Samp]) Frequeincy of sinusoid
f1=0.3;       % (Normalised [Hz/Samp]) Frequeincy of sinusoid
num_sig=500;  % Number of of realizations of random signal 

b=[1 1 1 1 1 1 1 1 1 1]; %filter coeficients (Low Pass MA filter)
a=0.1;

%_Noisy sinusoid_
for i=1:num_sig
    rng(i);      % Random generator seed     
    x(i,:) = a1*sin(f0*pi*2*(0:N-1)) + a2*sin(f1*pi*2*(0:N-1))+ sqrt(variance)*randn(1, N);
    x(i,:) = x(i,:)-mean(x(i,:));        % Zero mean
    acf_b=xcorr(x(i,:), 'biased');       % Biased ACF estimate
    acf_b=ifftshift(acf_b);
    xpb1(i,:)=fftshift(fft(acf_b));
end

%_filtered noisy sinusoid_
for i=1:num_sig
    rng(i);      % Random generator seed     
    x(i,:) = a1*sin(f0*pi*2*(0:N-1)) + a2*sin(f1*pi*2*(0:N-1))+ sqrt(variance)*randn(1, N);
    x(i,:) = filter(b, a, x(i,:));
    x(i,:) = x(i,:)-mean(x(i,:));        % Zero mean
    acf_b=xcorr(x(i,:), 'biased');       % Biased ACF estimate
    acf_b=ifftshift(acf_b);
    xpb2(i,:)=fftshift(fft(acf_b));
end

xpb_mean1(1,:)=mean(xpb1, 1); % Average of periodogram
xpb_var1(1,:)=var(xpb1, 1);   % variance of periodogram
xpb_mean2(1,:)=mean(xpb2, 1); % Average of periodogram
xpb_var2(1,:)=var(xpb2, 1);   % variance of periodogram

% plots of multiple realisiations of noisy sinusoid
subplot(3,2,1); hold on
plot(x_axis(2*N-1, 0.5), (xpb_mean1), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
plot(x_axis(2*N-1, 0.5), (xpb1),'color', [0.4010, 0.7450, 0.9830] )
plot(x_axis(2*N-1, 0.5), (xpb_mean1), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
title('PSD Estimates for Two Sinusoids in AWGN', 'FontWeight', 'normal');
ylabel('PSD');
xlim([0,0.5]); hold off
legend('Ensemble Average', 'PSD Estimates')

subplot(3,2,3); hold on
plot(x_axis(2*N-1, 0.5), 10*log10(xpb_mean1), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
plot(x_axis(2*N-1, 0.5), 10*log10(xpb1),'color', [0.4010, 0.7450, 0.9830] )
plot(x_axis(2*N-1, 0.5), 10*log10(xpb_mean1), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
title('PSD Estimates for Two Sinusoids in AWGN', 'FontWeight', 'normal');
ylabel('PSD [dB]');
ylim([-20, 35]); xlim([0,0.5]); hold off
legend('Ensemble Average', 'PSD Estimates')

% plot of variance of multiple realisiations of noisy sinusoid
subplot(3,2,5)
plot(x_axis(2*N-1, 0.5), sqrt(xpb_var1), 'lineWidth', 1.5, 'color', [0, 0.3470, 0.6410])
xlim([0,0.5])
title('Standard Deviation of the PSD Estimate', 'FontWeight', 'normal');
ylabel('Standard Deviation');


% plots of multiple realisiations of filtered noisy sinusoid
subplot(3,2,2); hold on
plot(x_axis(2*N-1, 0.5), (xpb_mean2), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
plot(x_axis(2*N-1, 0.5), (xpb2),'color', [0.4010, 0.7450, 0.9830] )
plot(x_axis(2*N-1, 0.5), (xpb_mean2), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
title('PSD Estimates for two Sinusoids in AWGN after LPF', 'FontWeight', 'normal');
ylabel('PSD');
xlim([0,0.5]); hold off
legend('Ensemble Average', 'PSD Estimates')

subplot(3,2,4); hold on
plot(x_axis(2*N-1, 0.5), 10*log10(xpb_mean2), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
plot(x_axis(2*N-1, 0.5), 10*log10(xpb2),'color', [0.4010, 0.7450, 0.9830] )
plot(x_axis(2*N-1, 0.5), 10*log10(xpb_mean2), 'lineWidth', 2.5, 'color', [0,0.3470, 0.6410])
title('PSD Estimates for two Sinusoids in AWGN after LPF', 'FontWeight', 'normal');
ylabel('PSD [dB]');
ylim([-25, 70]); xlim([0,0.5]); hold off
legend('Ensemble Average', 'PSD Estimates')

% plot of variance of multiple realisiations of filtered noisy sinusoid
subplot(3,2,6)
plot(x_axis(2*N-1, 0.5), sqrt(xpb_var2), 'lineWidth', 1.5, 'color', [0,0.3470, 0.6410])
title('Standard Deviation of the PSD Estimate', 'FontWeight', 'normal');
ylabel('Standard Deviation');
xlim([0,0.5])


% Common settings for all plots
for i=1:6
    subplot(3,2,i)
    xlabel('Normalised Frequency [cycles/sample]');
    box off; grid on
end

