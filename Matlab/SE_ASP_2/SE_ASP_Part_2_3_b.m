clc;clear;clf;
load EEG_Data_Assignment2.mat
N=length(POz);

%Pre-processing
POz=detrend(POz);

nfft=4096;                  % length of windowed segments
nstft=16384./4;                % (zero padding)
noverlap=round(nfft*0.75);  % Overlap percentage
win=hann(nfft);             % window

spectrogram(POz,win,noverlap,nstft,fs,'yaxis')
ylim([0,60])
title('Spectogram of POz EEG Data', 'fontWeight', 'normal')
xlabel('Time [Minutes]')
ylabel('Frequency [Hz]')
caxis([-140, -100])
