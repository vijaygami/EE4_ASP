clc;clear;
rng(1);

N=1000;                             % Length noise generated
a=[1, -2.76, 3.81, -2.65, 0.92];    % Model Coeficients
b=[1];

% AR process synthesis
wgn=randn(N,1); wgn = wgn-mean(wgn); wgn = wgn./sqrt(var(wgn));  % Zero mean, unit variance WGN
x=filter(b, a, wgn);                                             % AR process
x=x(501:end);                                                    % Remove transient part of IIR filter output


%p=[2:1:14];    % Model orders to plot for
p=[4, 6, 10];          % Model orders which yield best results
figure(1);clf; box off; grid on; hold on
for i=1:length(p)
    
    [a_est,o] = aryule(x,p(i));                         % Estimated model by yule Walker equations
    [h1,w]=freqz(sqrt(o), a_est, 'whole', 1024, 1);
    plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);    % Plotting estimates model estimates.

end

% Averaged periodogram welch method
[p_welch, w] = pwelch(x, hann(150), 0.9*150, 1024, 1);
plot(w, 10*log10(p_welch))

% True PSD from known model
[h1,w]=freqz(1, [1, -2.76, 3.81, -2.65, 0.92], 'whole', 1024, 1);  
plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1.5, 'color', 'black')
 
xlabel('Normalised Frequency [cycles/samp]');
ylabel('PSD [dB/cycles/samp]');
title('Spectrum of AR(4) Process of length 500. x(n) = 2.76x(n-1) - 3.81x(n-2) + 2.65x(n-3) - 0.92x(n-4) + w(n)', 'FontWeight', 'normal');
legend ('AR(4) Model Estimate', 'AR(6) Model Estimate','AR(10) Model Estimate', 'Averaged Periodogram (Welch)', 'True PSD', 'Location', 'best')
 
xlim([0,0.5]);ylim([-25, 50]); hold off;


%% pole zero plots of the estimate

figure(2);clf; box off; grid on; hold on
subplot(1,4,1)
zplane(b, a)
title({'Pole Zero plot of True ',
    '~~~~~~AR(4) Model'}, 'FontWeight', 'normal');

axis([-1.2 1.2 -1.2 1.2])
for i=1:length(p)
    
    [a_est,o] = aryule(x,p(i));                         % Estimated model by yule Walker equations
    subplot(1,4,i+1)
    zplane(sqrt(o), a_est)
    axis([-1.2 1.2 -1.2 1.2])
    title({'Pole Zero plot of',
        ['~AR(', num2str(p(i)), ') Estimate']}, 'FontWeight', 'normal');

end

    title({'Pole Zero plot of',
        ['AR(', num2str(p(i)), ') Estimate']}, 'FontWeight', 'normal');



