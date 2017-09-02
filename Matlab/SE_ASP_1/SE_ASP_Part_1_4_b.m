clc;clear; load EEG_Data_Assignment1.mat
N=length(POz);

%Pre-processing
POz=detrend(POz);

%% Blackman-Tukey with chebychev window with whole dataset,(very slow, 30sec, so commented out temprarily)

sidelobe=100;                              % Sidelobe attenuation in dB
M=round(N/10);                              % Window length for blackman tukey
L=2*N;                                     % FFT length, zero padd
x_cor = xcorr(POz', 'biased');              % Biased ACF

x_cor = x_cor(N-M:N+M);                              % Use 2M+1 values of the ACF
x_cor = x_cor.*chebwin(length(x_cor), sidelobe)';    % Chebychev window
x_cor = ifftshift(x_cor);
k=length(x_cor); k=(k+1)/2;
x_cor = [x_cor(1:k), zeros(1,L-2*k+1), x_cor(k+1:2*k-1)];  % Zero padd ACF while maintaining symmetry
pxx_bt=fftshift(fft(x_cor));                               % PSD estimate

figure(1); clf;
subplot(3,2,3); hold on
plot(x_axis(L,fs/2), 10*log10(pxx_bt))
hold off; xlim([0,60]); ylim([-110, -70]);

title('Blackman-Tukey Method', 'fontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('PSD [dB/Hz]')
grid on; box off;


%% Averaged periodograms

num_bins=10;            % Number of frequency bins per Hz
win_len=[1, 5, 10];     % window lengths in seconds

% Averaged over 1, 5, 10 seconds with hanning window
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    
    for j=0:(N/(win_len(i)*fs))-1
        est_seg(:, j+1)=pgm( hann(sl).*(POz( (j*sl+1) : sl*(j+1) ))  , L);    
    end
    
    pxx_ave_han(:,i)=mean(est_seg, 2);   % Average periodogram
    
end

% Averaged over 1, 5, 10 seconds with no window
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    
    for j=0:(N/(win_len(i)*fs))-1
        est_seg(:, j+1)=pgm( (POz( (j*sl+1) : sl*(j+1) ))  , L);    
    end
    
    pxx_ave(:,i)=mean(est_seg, 2);   % Average periodogram
    
end

%% Welch method


num_bins=10;            % Number of frequency bins per Hz
ovl_factor = 70;        % Percentage overlap
win_len=[1, 5, 10];     % window lengths in seconds

% Averaged over 1, 5, 10 seconds with hanning window
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((N-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(:, j+1)=pgm( hann(sl).*(POz( (j*step+1) : (j*step+sl)))  , L);    
    end
    
    pxx_welch(:,i)=mean(est_seg, 2);   % Average periodogram
    
end


%% Various Plots
subplot(3,2,1); hold on;
plot(x_axis(length(POz)*4, fs/2), 10*log10(pgm(POz, length(POz)*4)))

title('Periodogram', 'FontWeight', 'normal');
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
hold off; xlim([0,60]); ylim([-140, -70]); box off; grid on;

subplot(3,2,4)
plot(x_axis(L, fs/2), 10*log10(pxx_ave));hold on;
title('Averaged Periodogram', 'FontWeight', 'normal');
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
legend ('1 Second Window','5 Second Window','10 Second Window', 'location', 'best')
xlim([0,60]); ylim([-112, -65]); box off; grid on; hold off;

subplot(3,2,5)
plot(x_axis(L, fs/2), 10*log10(pxx_ave_han)); hold on
title('Averaged Periodogram, Hanning Window', 'FontWeight', 'normal');
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
legend ('1 Second Window','5 Second Window','10 Second Window', 'location', 'best')
xlim([0,60]); ylim([-112, -70]); box off; grid on; hold off

subplot(3,2,6)
plot(x_axis(L, fs/2), 10*log10(pxx_welch)); hold on
title('Averaged Periodogram, Welch Method with Hanning Window, $70\%$ overlap', 'FontWeight', 'normal');
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
legend ('1 Second Window','5 Second Window','10 Second Window', 'location', 'best')
xlim([0,60]); ylim([-112, -70]); box off; grid on; hold off


