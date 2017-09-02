clc;clear;rng(1)

% load data and pre-process
load('EEG_Data_Assignment2.mat')
POz=detrend(POz);
Cz=detrend(Cz);
N=length(POz);

% generate secondary noise input for ANC
o2=0.1;                      % Variance of noise          
std=sqrt(o2);                % Standard deviation of noise
wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));   % WGN input.
s=sin(2*pi*50/fs*[1:N])+wgn; % secondary noise input for ANC

% optimal periodogram parameters determined from 2.3.b
nfft=4096;                  % length of windowed segments
nstft=8192./2;              % (zero padding)
noverlap=round(nfft*0.75);  % Overlap percentage
win=hann(nfft);             % window

u_set=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05 0.1];   % Adapataion gains to try 
M_set=[5 10 20 50];       % Filter orders to try


 % plot of orginal POz
 figure(1); clf;
 subplot(2,4,1)
 spectrogram(POz,win,noverlap,nstft,fs,'yaxis')
 title('Unprocessed POz Spectrogram', 'fontWeight', 'normal')
 ylabel('Frequency [Hz]')
 xlabel('Time [Minutes]')
 ylim([0,120])
 caxis([-145, -100])
 
 
 % Varying adaptation gain
 x_est=zeros(length(u_set), N);
 for i =1:length(u_set)
     [x_est, w] = ANC_lms(POz, s, u_set(i), 10);
     subplot(2,4,i+1)
     spectrogram(x_est,win,noverlap,nstft,fs,'yaxis')
     ylim([0,120])
     title(['POz Spectrogram, M=10, $\mu$=', num2str(u_set(i))], 'fontWeight', 'normal')
     ylabel('Frequency [Hz]')
     xlabel('Time [Minutes]')
     caxis([-145, -100])
 end

% Remove colorbars
col=findall(gcf,'tag','Colorbar');
delete(col);
 
 
% Varying filter order
figure(2); clf;
for i =1:length(M_set)
    
    [x_est, w] = ANC_lms(POz, s, 0.001, M_set(i));
    subplot(2,length(M_set),i)
    spectrogram(x_est,win,noverlap,nstft,fs,'yaxis')
    ylim([0,120])
    title(['POz Spectrogram, M=', num2str(M_set(i)), ', $\mu$=0.001'], 'fontWeight', 'normal')
    ylabel('Frequency [Hz]')
    xlabel('Time [Minutes]')
    caxis([-145, -100])
end


% Remove colorbars
col=findall(gcf,'tag','Colorbar');
delete(col);



