clc;clear;
N=1500;
fs=1100;    % Sample frequency for FM sinewave
n_var=0.05; % AWCGN noise varaince

[y, fn] = fm_gen(fs, n_var);    % Generates FM signal

figure(1); clf; 
subplot(3,2,1)
plot(fn)
title('Instantaneous Frequency f(n)=d$\phi$/dn', 'fontWeight', 'normal')
ylabel('Instantaneous Frequency, f(n), [hz]')
xlabel('Time [Samples]')   
ylim([0,600]); 
box off; grid on;

subplot(2,2,2)
plot(real(y))
title('Real part of FM signal y(n)', 'fontWeight', 'normal')
ylabel('Real(y(n))'); 
xlabel('Time [Samples]')
ylim([-1.6,1.6]); 
box off; grid on;

subplot(2,2,4)
plot(imag(y))
title('Imaginary part of FM signal y(n)', 'fontWeight', 'normal')
ylabel('Imag(y(n))') 
xlabel('Time [Samples]')
ylim([-1.6,1.6]); 
box off; grid on;

% Find AR Model and plot PSD estimate
orders=[1, 3, 30];   % model orders to plot for
figure(2); %clf
for i=1:3
    
[a,o] = aryule(y,orders(i));
[h,w] = freqz(sqrt(o), a, 2048, fs, 'whole');

subplot(1,3,i)
plot(w, 10*log10(abs(h).^2));
title(['AR spectrum Estimate of y(n), p = ', num2str(orders(i))], 'fontWeight', 'normal')
xlabel('Frequency [Hz]')
ylabel('PSD Estimate [dB/Hz]')
grid on; box off;
xlim([0, fs])
end

mean(fn)    % mean instantaneous frequency

