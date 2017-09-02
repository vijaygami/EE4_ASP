clc;clear; figure(1);clf;

N_set=[1500, 10500];                      % Array of Signal lengths
a=[1, -2.76, 3.81, -2.65, 0.92];    % Model Coeficients
b=[1];

l=0;
for N = N_set
     l=l+1; rng(67);
    % AR process synthesis
    wgn=randn(N,1); wgn = wgn-mean(wgn); wgn = wgn./sqrt(var(wgn));  % Zero mean, unit variance WGN
    x=filter(b, a, wgn);                                             % AR process
    x=x(501:end);                                                    % Remove transient part of IIR filter output
    
    subplot(1, length(N_set), l);  box off; grid on; hold on;
    
    % True PSD from known model
    [h1,w]=freqz(1, [1, -2.76, 3.81, -2.65, 0.92], 'whole', 2048, 1);  
    plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1, 'color', 'black')
    
    p=[2,4, 200];          % Model orders which yield best results
    for i=1:length(p)
        [a_est,o] = aryule(x,p(i));                        % Estimated model by yule Walker equations
        [h1,w]=freqz(sqrt(o), a_est, 'whole', 2048, 1);
        plot(w, 10*log10(abs(h1).^2), 'lineWidth', 1);  % Plotting estimates model estimates.
    end

    xlabel('Normalised Frequency [cycles/sample]');
    ylabel('PSD [dB/cycles/sample]');
    title({['~~~~~~~~~~~~~~Spectrum of AR(4) Process of length ', num2str(N-500)], 
        'x(n) = 2.76x(n-1) - 3.81x(n-2) + 2.65x(n-3) - 0.92x(n-4) + w(n)'}, 'FontWeight', 'normal');
    legend ('True PSD', 'AR(2) Estimate', 'AR(4) Estimate','AR(200) Estimate', 'Location', 'best')
    xlim([0,0.5]);ylim([-25, 55]); hold off;
    
end

%% finding Average percentage error of coefficients and variance.
N_set=[1500, 10500, 100500, 1000500];
a_est=[];
for i=1:length (N_set)
    rng(1)
    wgn=randn(N_set(i),1); wgn = wgn-mean(wgn); wgn = wgn./sqrt(var(wgn));  % Zero mean, unit variance WGN
    x=filter(b, a, wgn);                                                    % AR process
    x=x(501:end);   
    [a_est(i,:), o(i,1)] = aryule(x,4);                                     % Estimated model by yule Walker equations

    % percentage errors in coefficients and variances
    mean_err(i) = mean(abs(a-a_est(i,:))./a)*100;
    noise_error(i) = abs(sqrt(o(i))-1)*100;
end





