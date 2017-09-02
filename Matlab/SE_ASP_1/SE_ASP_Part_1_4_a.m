% Sunspot Data
clc;clear;clf; load Sunspot.dat
years=Sunspot(:,1);
suns=Sunspot(:,2);

N=length(suns);
L=2^14;                  % Length of FFT (zero padd)

% No pre procesing
pxx_1=pgm(suns, L);

% Zero mean
pxx_2=pgm(suns-mean(suns), L);

% Detrend
pxx_3=pgm(detrend(suns), L);

% Detrend and Zero mean   % no need for this since detrend removes most of the mean.
%pxx_4=pgm(detrend(suns)-mean(detrend(suns)), L);

% log then zero mean
suns_log = Sunspot(:,2); 
suns_log(suns_log==0) = 1;  % set minimum sunspot count to be 1 rather than 0, since log(0) is infinite.
suns_log=log(suns_log);     % Log of data
suns_log=suns_log-mean(suns_log); % Zero mean
pxx_5=pgm(suns_log, L);

% Averaged Detrended periodogram
segments=8;     % Number of segments
sl=N/segments;  % Segment Length
est_seg=[];     % Will be filled with periodograms of sements
y=detrend(suns); % Preprocessing

for i=0:segments-1
    est_seg(:, i+1)=pgm( (y( (i*sl+1) : sl*(i+1) ))  , L);    
end
pxx_6=mean(est_seg, 2);   % Average periodogram

% Averaged Detrend periodogram (4 averages)
segments=4;     % Number of segments
sl=N/segments;  % Segment Length
est_seg=[];     % Will be filled with periodograms of sements


for i=0:segments-1
    est_seg(:, i+1)=pgm( (y( (i*sl+1) : sl*(i+1) ))  , L);    
end
pxx_7=mean(est_seg, 2);   % Average periodogram

% Averaged Detrend periodogram (16 averages)
segments=16;    % Number of segments
sl=N/segments;  % Segment Length
est_seg=[];     % Will be filled with periodograms of sements

for i=0:segments-1
    est_seg(:, i+1)=pgm( (y( (i*sl+1) : sl*(i+1) ))  , L);    
end
pxx_8=mean(est_seg, 2);   % Average periodogram


%% Various plots of pre-processed sun data periodograms

% Time series and log of time series.
subplot (2, 2, 1)
hold on
plot(years, (suns-mean(suns))./sqrt(var(suns-mean(suns))));
plot(years, (suns_log)./sqrt(var(suns_log)));
hold off
ylim([-3.5, 5.5])
title('Sunspot Time Series', 'FontWeight', 'normal');
xlabel('Year');
ylabel('Normalised Count');
legend('Normalised Sunspot','Normalised Log of Sunspot', 'location', 'best')
box off; grid on; 

% Different Preprocesisng methods
subplot (2, 2, 2)
hold on
plot(x_axis(L, 0.5), 10*log10(pxx_1))
plot(x_axis(L, 0.5), 10*log10(pxx_2))
plot(x_axis(L, 0.5), 10*log10(pxx_3))
hold off; box off; grid on; 
xlim([0, 0.4]); ylim([-30,60])

title('Periodogram', 'FontWeight', 'normal');
xlabel('Frequency [cycles/year]');
ylabel('PSD [dB/cycles/year]');
legend('No Preprocessing', 'Zero Mean','Detrend', 'location', 'best')



% Periodogram of log of data
subplot (2, 2, 3);hold on
plot(x_axis(L, 0.5), 10*log10(pxx_5))
hold off; box off; grid on; 
title('Periodogram of Zero Mean Log of Data', 'FontWeight', 'normal');
xlabel('Frequency [cycles/year]');
ylabel('PSD [dB/cycles/year]');
legend('Log and Zero Mean', 'location', 'best')
ylim([-40,30]); xlim([0, 0.4]);

% plots of averaged periodograms
subplot (2, 2, 4)
hold on
plot(x_axis(L, 0.5), 10*log10(pxx_6))
plot(x_axis(L, 0.5), 10*log10(pxx_7))
plot(x_axis(L, 0.5), 10*log10(pxx_8))
hold off; box off; grid on; 
xlim([0, 0.4]); ylim([10,45]);

title('Averaged Periodogram', 'FontWeight', 'normal');
xlabel('Frequency [cycles/year]');
ylabel('PSD [dB/cycles/year]');
legend('4 Averages', '8 averages','16 Averages', 'location', 'best')



