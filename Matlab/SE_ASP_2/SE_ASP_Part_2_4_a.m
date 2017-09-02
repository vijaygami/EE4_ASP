%% Load Data
clc;clear;
fs=4;


%Vijay Gami Data
load vgtr1.mat
load vgtr2.mat
load vgtr3.mat

%{
%Will Smith Data set 1
load ws1tr1.mat
load ws1tr2.mat

%Will Smith Data set 2
load ws2tr1.mat
load ws2tr2.mat
load ws2tr3.mat
%}

%% Preprocessing
RRI_trial_1=detrend(RRI_trial_1);
RRI_trial_2=detrend(RRI_trial_2);
RRI_trial_3=detrend(RRI_trial_3);

%% Periodograms of whole data. Hanning window used.

figure(1)
subplot(1,3,1)
plot(x_axis(length(RRI_trial_1), fs/2), 10*log10(pgm(RRI_trial_1.*hann(length(RRI_trial_1))', length(RRI_trial_1))))
title('Periodogram for RRI Trial 1', 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency (dB/Hz)')
xlim([0,fs/2]);ylim([-90,-5]); box off; grid on

subplot(1,3,2)
plot(x_axis(length(RRI_trial_2), fs/2), 10*log10(pgm(RRI_trial_2.*hann(length(RRI_trial_2))', length(RRI_trial_2))))
title('Periodogram for RRI Trial 2', 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency (dB/Hz)')
xlim([0,fs/2]);ylim([-90,-5]); box off; grid on

subplot(1,3,3)
plot(x_axis(length(RRI_trial_3), fs/2), 10*log10(pgm(RRI_trial_3.*hann(length(RRI_trial_3))', length(RRI_trial_3))))
title('Periodogram for RRI Trial 3', 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency (dB/Hz)')
xlim([0,fs/2]);ylim([-90,0]); box off; grid on

%% Averaged periodograms for the three trials

num_bins=1000;                % Number of frequency bins per Hz
ovl_factor = 0;               % Percentage overlap
win_len=[12.5, 25, 37.5];     % window lengths in seconds. (corresponds to 50, 100, 150 samples)

% Compute averaged periodogram for Trial 1
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_1)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_1( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr1(i,:)=mean(est_seg, 1);   % Average periodogram

end

% Compute averaged periodogram for Trial 2
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_2)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_2( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr2(i,:)=mean(est_seg, 1);   % Average periodogram

end

% Compute averaged periodogram for Trial 3
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_3)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_3( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr3(i,:)=mean(est_seg, 1);   % Average periodogram

end

% averaged plots
figure(3); clf;
subplot(3, 2, 1);
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr1));
title(['Averaged Periodogram for RRI Trial 1, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;
  
subplot(3, 2, 3)
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr2));
title(['Averaged Periodogram for RRI Trial 2, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;

subplot(3, 2, 5)
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr3));
title(['Averaged Periodogram for RRI Trial 3, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;

%% same plots but with now with welch method 80% overlap

num_bins=1000;                % Number of frequency bins per Hz
ovl_factor = 80;               % Percentage overlap
win_len=[12.5, 25, 37.5];     % window lengths in seconds. (corresponds to 50, 100, 150 samples)

% Compute averaged periodogram for Trial 1
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_1)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_1( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr1(i,:)=mean(est_seg, 1);   % Average periodogram

end

% Compute averaged periodogram for Trial 2
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_2)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_2( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr2(i,:)=mean(est_seg, 1);   % Average periodogram

end

% Compute averaged periodogram for Trial 3
for i=1:length(win_len)
    
    sl=fs*win_len(i);  % Segment Length
    L=fs*num_bins;     % Length of FFT in order to have 'num_bins' per Hz
    est_seg=[];        % Will be filled with periodograms of sements
    step=(100-ovl_factor)*sl/100; % Step length
    
    segs=floor((length(RRI_trial_3)-sl)/step); % number of segments
    
    for j=0:segs-1
        est_seg(j+1, :) = pgm( hann(sl)'.*(RRI_trial_3( (j*step+1) : (j*step+sl))) , L );    
    end
    
    pxx_welch_tr3(i,:)=mean(est_seg, 1);   % Average periodogram

end

% averaged plots
figure(3); 
subplot(3, 2, 2);
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr1));
title(['Averaged Periodogram for RRI Trial 1, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;
  
subplot(3, 2, 4)
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr2));
title(['Averaged Periodogram for RRI Trial 2, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;

subplot(3, 2, 6)
plot(x_axis(L, fs/2), 10*log10(pxx_welch_tr3));
title(['Averaged Periodogram for RRI Trial 3, Overlap = ', num2str(ovl_factor),'$\%$'], 'FontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')
legend('50 Samples', '100 Samples', '150 Samples')
xlim([0, 2]); box off; grid on;

