clc;clear;

N=256;  % Length of signal
L=16384; % length of FFT

wb=(1/N)*[[N:-1:1],zeros(1,(L-(2*N-1))),[1:1:N-1]];  % Zero padded Bartlett window. Length N, padd to L
wbf=fftshift(abs(fft(wb, L)));

% Plot spectra of window
clf
hold on
plot(x_axis(L, 0.5), 20*log10(wbf));
xlim([0, 0.06])

plot([4, 4]./256,[-100, 20], 'LineWidth', 1.2, 'LineStyle', '--' )
plot([12, 12]./256,[-100, 0], 'LineWidth', 1.2, 'LineStyle', '--' )

% Adding the freqiencies coresponding to the two alphas
text(4/256,25,['$\alpha$ = 4'], 'FontSize', 12 )
text(12/256, 5,['$\alpha$ = 12'],'FontSize', 12 )

title(['Spectra of Bartlett Window of Length 256'], 'FontWeight', 'normal');
xlabel('Normalised Frequency [x2$\pi$ rad/samp]')
ylabel('Amplitude [dB]');
grid on
hold off