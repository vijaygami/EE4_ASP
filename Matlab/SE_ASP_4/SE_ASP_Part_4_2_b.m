clc;clear;
L=1500;     % length of FM signal    
fs=1100;    % Sample frequency for FM sinewave
n_var=0.05; % AWCGN noise varaince


y = fm_gen(fs, n_var);
[~, a, ~ ] = CLMS_AR(y, 0.075, 1);    % Run complex-valued LMS algorithm to estimate AR coefficient ^a1(n)


for n = 1:L
    [h ,w]= freqz(1 , [1; -conj(a(n))], 2048, fs); 	% Compute power spectrum for each coef
    H(:, n) = abs(h).^2;                            % Store it in a matrix
end

% Remove outliers in the matrix H
medianH = 50*median(median(H));
H(H > medianH) = medianH;
% Plot time-frequency diagram

figure(1); clf;
subplot(1,2,1)
surf(1:length(y), w./1000, H.^2, 'LineStyle', 'none')        % plot the adaptive LMS-AR based spectrum estimate
view(2)
xlabel('Iteration [n]', 'interpreter', 'latex')
ylabel('Frequency [Khz]', 'interpreter', 'latex')
title('CLMS AR(1) Time Frequency Spectrum of FM signal y(n)', 'fontWeight', 'normal')
ylim([0 0.55])

subplot(1,2,2)                                      % plot Spectrogramto compare
spectrogram(y, hann(128), 120, 1024, fs, 'yaxis');
view(2)
title('Spectrogram of FM signal y(n)', 'fontWeight', 'normal')
xlabel('Time [Seconds]', 'interpreter', 'latex')
ylabel('Frequency [KHz]', 'interpreter', 'latex')
ylim([0, 0.55])
caxis([-70, -20])


