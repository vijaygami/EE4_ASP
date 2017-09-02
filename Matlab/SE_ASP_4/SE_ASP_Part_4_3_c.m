clc;clear;
fs=1100;     % Sample frequency for FM sinewave
N=1100;      % Number of weights for DFT-CLMS   
n_var=0.05;  % AWCGN noise varaince

[y, fn] = fm_gen(fs, n_var);  % generate FM signal of length 1500

%% Apply DFT-CLMS

[ ~, w , ~] = DFT_CLMS(y, N, 0, 1);            % non-leaky
[ ~, wleak , ~] = DFT_CLMS(y, N, 0.05, 1);     % leaky

figure(1); clf;
subplot(1,2,1)
surf(1:length(y), ([0:N-1].*fs/N)./1000, abs(w), 'LineStyle', 'none')
view(2); ylim([0,0.55])
xlabel('Time [Samples]','interpreter', 'latex')
ylabel('Frequency [KHz]','interpreter', 'latex')
title('DFT-CLMS Applied to FM signal', 'fontWeight', 'normal')

subplot(1,2,2)
surf(1:length(y), ([0:N-1].*fs/N)./1000, abs(wleak).^2, 'LineStyle', 'none')
view(2); ylim([0,0.55])
xlabel('Time [Samples]', 'interpreter', 'latex')
ylabel('Frequency [KHz]', 'interpreter', 'latex')
title('Leaky DFT-CLMS Applied to FM signal, $\gamma$ = 0.05', 'fontWeight', 'normal')

%% Compare FFT and DFT-CLMS at time instant N

figure(2);clf
subplot(1,2,1); 
plot(abs(w(:,N+1)));
xlabel('Frequency [Hz]', 'interpreter', 'latex')
ylabel('$|w(N)|$', 'interpreter', 'latex')
title('DFT-CLMS Filter Weight at Time Index N, $\mu$=1', 'fontWeight', 'normal')
box off; grid on; xlim([0,1100])

subplot(1,2,2)
fft_weights=((1/N)*fft(y, N))';
plot(abs(fft_weights));
xlabel('Frequency [Hz]', 'interpreter', 'latex')
ylabel('(1/N)*$|FFT(y)|$', 'interpreter', 'latex')
title('FFT of length N', 'fontWeight', 'normal')
box off; grid on; xlim([0,1100])

MSWE=mean(abs(fft_weights-w(:,N+1)).^2)    % MSE in weights

%%
figure(3)
subplot(1,2,1)
[hAx,hLine1,hLine2] = plotyy([1:1500], abs(w((100*fs/N+1), :)), [1:1500], abs(w((350*fs/N+1), :)));    %100 Hz weight
legend('100 Hz Weight','350 Hz Weight')
title('Evolution of Weights of DFT-CLMS', 'fontWeight', 'normal')
xlabel('Time [Samples]', 'interpreter', 'latex')
ylabel(hAx(1), 'Weight Magnitude', 'color', 'black', 'interpreter', 'latex')
box off; grid on
set(gca,'xtick',[0:250:1500])
set(gca,'ytick',[0:0.1:1])
ylim([0,0.75])

subplot(1,2,2)
w_der=abs(w);
w_der=diff(w_der, 1, 2);
surf(1:(length(y)-1), ([0:N-1].*fs/N)./1000, w_der, 'LineStyle', 'none')
view(2); ylim([0,0.55])
xlabel('Time [Samples]', 'interpreter', 'latex')
ylabel('Frequency [KHz]', 'interpreter', 'latex')
title('Time Derivitive the DFT-CLMS Weigths', 'fontWeight', 'normal')


