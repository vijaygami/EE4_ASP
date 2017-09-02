%% load data and pre-process
clc;clear;
load('EEG_Data_Assignment2.mat')
POz=detrend(POz);
Cz=detrend(Cz);

a=1500;     % Starting point of data to use
L=1200;     % length of data to use
POz=POz(a:a+L-1);
N=5000;


%% Apply DFT_CLMS to POz and plot

figure(1); clf; subplot(1,2,1)
[ ~, w, ~ ] = DFT_CLMS(POz, N, 0, 1);
surf(1:L, [0:N-1].*fs/N, 10*log10(abs(w).^2), 'LineStyle', 'none')
view(2)
xlabel('Iteration [n]')
ylabel('Frequency [Hz]')
title('DFT-CLMS Applied to Real World Poz Data', 'fontWeight', 'normal')
ylabel(colorbar, 'Power/Frequency [dB/Hz]')
ylim([0,60])
xlim([0, L])

subplot(1,2,2)
surf(1:L, [0:N-1].*fs/N, 10*log10(abs(w).^2), 'LineStyle', 'none')
view(2)
xlabel('Iteration [n]')
ylabel('Frequency [Hz]')
title('DFT-CLMS Applied to Real World Poz Data', 'fontWeight', 'normal')
ylim([0,60])
xlim([0, L])
caxis([-137, -110]) 
ylabel(colorbar, 'Power/Frequency [dB/Hz]')
