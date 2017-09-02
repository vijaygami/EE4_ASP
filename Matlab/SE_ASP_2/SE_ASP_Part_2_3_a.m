clc;clear;clf
t = -5:0.001:5;
fs=2000;

x = chirp(t,100,1,200,'quadratic');


%% investigate effects of zero padding (FFT length)
% 64 window length, no overlap

figure(1)
nfft=64;               % length of windowed segments
nstft=64;             % (zero padding)
noverlap=round(nfft*0); % percentage overlap
subplot(4,4,1);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=64, FFT len=64', 'fontWeight', 'normal')

nfft=64;               % len of windowed segments
nstft=128;             % (zero padding)
noverlap=round(nfft*.0); % percentage overlap
subplot(4,4,2);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=64, FFT len=128', 'fontWeight', 'normal')

nfft=64;               % length of windowed segments
nstft=256;             % (zero padding)
noverlap=round(nfft*.0); % percentage overlap
subplot(4,4,3);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=64, FFT len=256', 'fontWeight', 'normal')

nfft=64;               % length of windowed segments
nstft=512;             % (zero padding)
noverlap=round(nfft*.0); % percentage overlap
subplot(4,4,4);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=64, FFT len=512', 'fontWeight', 'normal')



%% investigate effects of window length
% no window, no overlap

nfft=64;                % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*0); % percentage overlap
subplot(4,4,5);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=64, FFT len=64', 'fontWeight', 'normal')


nfft=128;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*0); % percentage overlap
subplot(4,4,6);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=128, FFT len=128', 'fontWeight', 'normal')


nfft=256;               % len of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*0); % percentage overlap
subplot(4,4,7);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=256, FFT len=256', 'fontWeight', 'normal')


nfft=512;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*0); % percentage overlap
subplot(4,4,8);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Window len=512, FFT len=512', 'fontWeight', 'normal')


%% investigate effects of percentage overlap

nfft=256;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*.2); % percentage overlap
subplot(4,4,9);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Percentage Overlap=20$\%$', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*.6); % percentage overlap
subplot(4,4,10);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Percentage Overlap=60$\%$', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*.8); % percentage overlap
subplot(4,4,11);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Percentage Overlap=80$\%$', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=nfft;             % (zero padding)
noverlap=round(nfft*.9); % percentage overlap
subplot(4,4,12);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Percentage Overlap=90$\%$', 'fontWeight', 'normal')



%% investigate effects of different windows

nfft=256;               % length of windowed segments
nstft=2048;             % (zero padding)
noverlap=round(nfft*.9); % percentage overlap
subplot(4,4,13);
spectrogram(x,rectwin(nfft),noverlap,nstft,fs,'yaxis')
title('Rectangular Window', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=2048;             % (zero padding)
noverlap=round(nfft*.9); % percentage overlap
subplot(4,4,14);
spectrogram(x,hamming(nfft),noverlap,nstft,fs,'yaxis')
title('Hamming Window', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=2048;             % (zero padding)
noverlap=round(nfft*.9); % percentage overlap
subplot(4,4,15);
spectrogram(x,hann(nfft),noverlap,nstft,fs,'yaxis')
title('Hanning Window', 'fontWeight', 'normal')

nfft=256;               % length of windowed segments
nstft=2048;             % (zero padding)
noverlap=round(nfft*.9); % percentage overlap
subplot(4,4,16);
hold on
spectrogram(x,chebwin(nfft, 150),noverlap,nstft,fs,'yaxis')
title('Chebyshev Window, 150dB Att.', 'fontWeight', 'normal')
% Ploting the actual instantaneous frequency of the chirp
plot(t+2, 100+100*t.^2, 'lineWidth', 1, 'color', 'red')
xlim([0.1,3.9])

% hide all colorbars
ch=findall(gcf,'tag','Colorbar');
delete(ch);

for i=1:16
    subplot(4,4,i)
    ylabel('Frequency [Hz]')
    xlabel('Time [Seconds]')
end

