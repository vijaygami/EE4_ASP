clc;clear;clf;

rng(5);      % Random generator seed     
N=1000;      % Length of signal
variance=1;  % Variance of AWGN
a1=2;        % Amplitude of sinusoid
f0=0.2;      % Normalised [Hz/Samp]) Frequency of sinusoid

b=[1 1 1 1 1 1 1 1 1 1]; % Filter coeficients (MA LPF filter chosen)
a=[0.1];

%% _Noisy Sinusoid_
x = a1*sin(f0*pi*2*(0:N-1)) + sqrt(variance)*randn(1, N); x=x-mean(x);  % zero mean

[acf_b, acf_ub]=computeACF(x); % Computes biased and unbiased ACF estimate

acf_b=ifftshift(acf_b);     % biased ACF
acf_ub=ifftshift(acf_ub);   % unbiased ACF

px_b=fftshift(fft(acf_b));
px_ub=fftshift(fft(acf_ub));

% ACF plots
subplot(2,3,1); box off; grid on; hold on
plot(acf_ub);
plot(acf_b);
legend('Unbiased', 'Biased', 'location', 'best');
title('ACF estimate of Noisy Sinusoid', 'FontWeight', 'normal');
ylabel('ACF')
xlabel('Correlation Lag');
hold off
xlim([0, N]); ylim([-3.5,4.5])

% Correlogram plots
subplot(2,3,4); box off; grid on; hold on
plot(x_axis(length(px_ub), 0.5), px_ub);
plot(x_axis(length(px_b), 0.5), px_b);
legend('Unbiased', 'Biased', 'location', 'best')
hold off
title('Correlogram of Noisy Sinusoid', 'FontWeight', 'normal');
xlabel('Normalised Frequency [cycles/sample]')
ylabel('PSD Estimate');
xlim([0.175, 0.225]);

%% _WGN_
x = randn(1, N);
x=x-mean(x);

[acf_b, acf_ub]=computeACF(x); %computes biased and unbiased ACF estimate

acf_b=ifftshift(acf_b);
acf_ub=ifftshift(acf_ub);

px_b=fftshift(fft(acf_b));
px_ub=fftshift(fft(acf_ub));

% ACF plots
subplot(2,3,2); box off; grid on; hold on
plot(acf_ub);
plot(acf_b);
legend('Unbiased', 'Biased', 'location', 'best')
title('ACF estimate of WGN', 'FontWeight', 'normal');
ylabel('ACF')
xlabel('Correlation Lag');
xlim([0, N]); hold off;

% Correlogram plots
subplot(2,3,5); box off; grid on; hold on
plot(x_axis(length(px_ub), 0.5), (px_ub));
plot(x_axis(length(px_b), 0.5), (px_b));
legend('Unbiased', 'Biased', 'location', 'best')
title('Correlogram of WGN', 'FontWeight', 'normal');
xlabel('Normalised Frequency [cycles/sample]')
ylabel('PSD Estimate');
xlim([0, 0.5]); ylim([-30,30]); hold off


%% _Filtered WGN_
x = randn(1, N);
x=x-mean(x);
x=filter(b,a,x);   

[acf_b, acf_ub]=computeACF(x); %computes biased and unbiased ACF estimate

acf_b=ifftshift(acf_b);
acf_ub=ifftshift(acf_ub);

px_b=fftshift(fft(acf_b));
px_ub=fftshift(fft(acf_ub));

% ACF plots
subplot(2,3,3); box off; grid on; hold on
plot(acf_ub);
plot(acf_b);
legend('Unbiased', 'Biased', 'location', 'best')
hold off
title('ACF Estimate of Filtered WGN', 'FontWeight', 'normal');
ylabel('ACF')
xlabel('Correlation Lag');
xlim([-1, N])

% Correlogram plots
subplot(2,3,6); box off; grid on; hold on
plot(x_axis(length(px_ub), 0.5), (px_ub));
plot(x_axis(length(px_b), 0.5), (px_b));
legend('Unbiased', 'Biased', 'location', 'best')
title('Correlogram of Filtered WGN', 'FontWeight', 'normal');
xlabel('Normalised Frequency [cycles/sample]')
ylabel('PSD Estimate');
hold off
xlim([0, 0.2])
%ylim([-40,40])

