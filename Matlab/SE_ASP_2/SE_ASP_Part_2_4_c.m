%% Load Data
clc;clear;
fs=4;

%Vijay Gami Data
load vgtr1.mat
load vgtr2.mat
load vgtr3.mat

%% Preprocessing
RRI_trial_1=detrend(RRI_trial_1);
RRI_trial_2=detrend(RRI_trial_2);
RRI_trial_3=detrend(RRI_trial_3);

%% PACF to determine model orders.
figure(1); clf;
subplot(1,3,1)
parcorr(RRI_trial_1, 40)

subplot(1,3,2)
parcorr(RRI_trial_2, 40)

subplot(1,3,3)
parcorr(RRI_trial_3, 40)

for i=1:3
    subplot(1,3,i)
    title(['Partial ACF for ECG data, Trial ', num2str(i)], 'FontWeight', 'normal')
    xlabel('Correlation lag')
    grid on
    box off
    legend('PACF', '$95\%$ Confidence Interval')
end

%% MDL AIC to determine optimal model orders.
% Code adapted from my (Vijay Gami) ASP CW year 3. 

or=40;       %max order
p=[1:or];    %Model orders

PACF_1=[]; %will store store partial ACF for trial 1
for i=1:or,
    [a,e1(i)]=aryule(RRI_trial_1,i);
end

PACF_2=[]; %will store store partial ACF for trial 2
for i=1:or,
    [a,e2(i)]=aryule(RRI_trial_2,i);
end

PACF_3=[]; %will store store partial ACF for trial 3
for i=1:or,
    [a,e3(i)]=aryule(RRI_trial_3,i);
end

figure(2); clf;
subplot(1,2,1); % MDL
MDL1=log10(e1)+p.*(log10(length(RRI_trial_1))/length(RRI_trial_1)); % Trial 1
MDL2=log10(e2)+p.*(log10(length(RRI_trial_2))/length(RRI_trial_2)); % Trial 2
MDL3=log10(e3)+p.*(log10(length(RRI_trial_3))/length(RRI_trial_3)); % Trial 3
plot(p,MDL1, p, MDL2, p, MDL3)
legend('MDL Trial 1', 'MDL Trial 2', 'MDL Trial 3')
title('MDL for the Three Trials', 'FontWeight', 'normal')
xlabel('Model Order P')
ylabel('MDL')
box off, grid on;

subplot(1,2,2); % AIC
AIC1=log10(e1)+p.*(2/length(RRI_trial_1)); % Trial 1
AIC2=log10(e2)+p.*(2/length(RRI_trial_2)); % Trial 2
AIC3=log10(e3)+p.*(2/length(RRI_trial_3)); % Trial 3
plot(p, AIC1, p, AIC2, p, AIC3)
legend('AIC Trial 1', 'AIC Trial 2', 'AIC Trial 3')
title('AIC for the Three Trials', 'FontWeight', 'normal')
xlabel('Model Order P')
ylabel('AIC')
box off, grid on;


%% plotting AR Models
figure(3); clf;
winlen=128; % window length for pwelch method)
nfft=2048;  % FFT length (zero padd)


% Trial 1
subplot(2,3,1); hold on
[a,or] = aryule(RRI_trial_1,2);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_1,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (Minimum Order) Model for Trial 1', 'FontWeight', 'normal')
legend('AR(2) Model', 'Welch Periodogram')
hold off

subplot(2,3,4); hold on
[a,or] = aryule(RRI_trial_1,22);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_1,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (High Order) Model for Trial 1', 'FontWeight', 'normal')
legend('AR(22) Model', 'Welch Periodogram')

hold off

% Trial 2
subplot(2,3,2); hold on
[a,or] = aryule(RRI_trial_2,9);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_2,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (Minimum Order) Model for Trial 2', 'FontWeight', 'normal')
legend('AR(9) Model', 'Welch Periodogram')
hold off

subplot(2,3,5); hold on
[a,or] = aryule(RRI_trial_2,19);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_2,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (High Order) Model for Trial 2', 'FontWeight', 'normal')
legend('AR(19) Model', 'Welch Periodogram')

hold off

% Trial 3
subplot(2,3,3); hold on
[a,or] = aryule(RRI_trial_3,2);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_3,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (Minimum Order) Model for Trial 3', 'FontWeight', 'normal')
legend('AR(2) Model', 'Welch Periodogram')

hold off

subplot(2,3,6); hold on
[a,or] = aryule(RRI_trial_3,34);                                 % AR model
[h1,w]=freqz(sqrt(or), a, nfft, fs);
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);
[pxx w]=pwelch(RRI_trial_3,hann(winlen),winlen*0.8,nfft, fs);   % welch periodogram to compare
plot(w, 10*log10(pxx));
title('AR (High Order) Model for Trial 3', 'FontWeight', 'normal')
legend('AR(34) Model', 'Welch Periodogram')

hold off

% Common plot settings
for i=1:6
    subplot(2,3,i)
    xlim([0,2])
    box off; grid on
    xlabel('Frequency [Hz]')
    ylabel('PSD [dB/Hz]')
end


