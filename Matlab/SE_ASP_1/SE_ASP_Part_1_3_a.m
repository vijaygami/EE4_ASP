clc;clear;clf;

% The set of values of windw lengths (N) to compute and plot sidelobe peaks and 3dB bandwidth for.
% Note total window length is 2N.
set=[[6:2:32], 64, 128, 256, 512, 1024, 2048, 4096]; 

oversamp = 2048; % Ratio of FFT length to N (zero pad)
Num_lobes = 3;   % Number of sidelobe peaks to plot

% Empty vectors will be filled with computed values of sidelobe peaks and 3dB BW.
peaks_set=[];   
BW_set=[];
error=[];       % Stores error of computed 3dB point so it can be monitered for accuracy reasons
for N = set

    L=N*oversamp+1; %FFT length
    window=(1/N)*[[N:-1:1],zeros(1,(L-(2*N-1))), [1:1:N-1]];  % Zero padded bartlett window of length L
    f=abs(fft(window));

    [~, index] = min(abs(f - N/(sqrt(2))));    % Find the closest 3dB point
    norm_freq=(index-1)*pi*2/L;                % 3dB bandwidth in normalised frequency [rad/sample]
    
    error = [error, abs(1/sqrt(2) - (f(index)/N))]; % Error between approx 3dB point
    
    % Find attenuation of sidelobes 
    peaks = findpeaks(f);
    peaks=20*log10(peaks(1:Num_lobes)/N)';
    
    % Append to array the computed values
    BW_set=[BW_set, norm_freq];                 
    peaks_set=[peaks_set, peaks];
    
end

% Plot 3dB points vs N 
figure(1)
subplot(1, 3, 1)
semilogx(set,BW_set,'LineWidth', 1.2); %x-axis is log scale
xlabel('N')
ylabel('3dB BW [rad/sample]');
title('Mainlobe 3dB Width', 'FontWeight', 'normal');
grid on
box off

% Plot 3dB points vs N (log-log)
subplot(1, 3, 2)
plot(log10(set),log10(BW_set),'LineWidth', 1.2)
xlabel('log10(N)')
ylabel('log10(3dB BW) [log(rad/sample)]');
title('Mainlobe 3dB Width (log10)', 'FontWeight', 'normal');
grid on
box off

% Linear line of best fit
p = polyfit(log10(set),log10(BW_set),1);
text(0.7,-.3,['log10 (BW) = ', num2str(p(1,1)), '*log10 (N) + ', num2str(p(1,2))] )

% Plot sidelobe ampltudes vs N
subplot(1, 3, 3)
semilogx(set,peaks_set,'LineWidth', 1.2)
xlabel('N')
ylabel('Sidelobe gain [dB]');
title('Sidelobe Gain [dB]', 'FontWeight', 'normal');
legend('Sidelobe 1','Sidelobe 2','Sidelobe 3', 'Location','east');
grid on
box off




%_Verifying results through some (normalised) plots of the window__________________________________________________
figure (2)
set=[128, 256, 512]; 

oversamp = 2048; % Ratio of FFT length to N (zero pad)

i=0;
for N = set
    i=i+1;
    
    L=N*oversamp+1; %FFT length
    window=(1/N)*[[N:-1:1],zeros(1,(L-(2*N-1))), [1:1:N-1]];  % Zero padded bartlett window of length L
    f=abs(fft(window))./N;  % Normalised Amplitude for convineince

    [~, index] = min(abs(f - 1/(sqrt(2))));    % Find the closest 3dB point
    norm_freq=(index-1)*pi*2/L;                % 3dB bandwidth in normalised frequency [rad/sample]
    
    
    % Find attenuation of sidelobes 
    [peaks, locs] = findpeaks(f);
    peak=20*log10(peaks(1,1))';
    loc=(locs(1,1)-1)*pi*2/L; 
    
    % Log plots

    subplot(2, 3, i)
    hold on
    plot(x_axis(L, pi),20*log10(fftshift(f)), 'LineWidth', 1.2);
    plot([0, norm_freq], [-3, -3], 'color', [0.5 0 0])
    plot([norm_freq, norm_freq], [-100, -3], 'color', [0.5 0 0])
    plot([0, loc], [peak, peak], 'color', [0 .5 0])
    plot([loc, loc], [-100, peak], 'color', [0 .5 0])
    
    plot(norm_freq,-3, 'o', 'color', [0.5 0 0])    
    plot(loc,peak, 'o', 'color', [0 0.5 0])
    
    text(norm_freq,-1,['(', num2str(norm_freq), ', ', num2str(-3), ')'] )
    text(loc-0.005,peak+2,['(', num2str(loc), ', ', num2str(peak), ')'] )
    hold off
    
    xlim([0, 0.08])
    ylim([-40, 1])
    xlabel({'Normalised Frequency [rad/samp]'});
    ylabel('Amplitude [dB]');
    title(['Bartlett Window 3dB BW and Sidelobe, N = ',num2str(N)], 'FontWeight', 'normal');
    grid on

 
    % linear plots
    subplot(2, 3, i+3)
    hold on
    plot(x_axis(L, pi),fftshift(f), 'LineWidth', 1.2);
    plot([0, norm_freq], [1./sqrt(2), 1./sqrt(2)], 'color', [0.5 0 0])
    plot([norm_freq, norm_freq], [0, 1./sqrt(2)], 'color', [0.5 0 0])
    plot([0, loc], [10^(peak/20), 10^(peak/20)], 'color', [0 .5 0])
    plot([loc, loc], [0, 10^(peak/20)], 'color', [0 .5 0])
    
    plot(norm_freq,1./sqrt(2), 'o', 'color', [0.5 0 0])    
    plot(loc,10^(peak/20), 'o', 'color', [0 0.5 0])
    
    text(norm_freq,0.75,['(', num2str(norm_freq), ', ', num2str(0.7071), ')'] )
    text(loc,10^(peak/20)+0.05,['(', num2str(loc), ', ', num2str(10^(peak/20)), ')'] )
    hold off
    grid on 
    
    xlim([0, 0.08])
    ylim([0, 1]);
    xlabel({'Normalised Frequency [rad/samp]'});
    ylabel('Amplitude');
    title(['Bartlett Window 3dB BW and Sidelobe, N = ',num2str(N)], 'FontWeight', 'normal');
end



