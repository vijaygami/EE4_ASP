clc;clear;

N=1000;             % Number of input samples
wo=0.01*pi;         % Frequency of input signal
o2=1;               % Variance of noise
std=sqrt(o2);       % Standard deviation of noise

M=[1:1:25];         % Filter  order;
u=0.01;             % adaptation gain of ANC LMS
num_ens=100;        % number of independnat trials

x=sin(wo.*[1:N]);   % Signal generation

x_est=zeros(num_ens,N); % pre-allocate memory for speed
MSPE=zeros(num_ens,1);
MSPE_set_ANC=zeros(1,length(M));
for j=1:length(M)
    for i = 1:num_ens
        rng(i)
        wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));    % WGN input.
        eta=filter([1 0 0.5],[1], wgn);                                    % Noise Generation;
        s=x+eta;                                                           % Noisy signal

        [x_est(i,:), w] = ANC_lms(s, wgn, u, M(j));
        MSPE(i,1) = mean((x_est(i,:) - x).^2);
    end

    MSPE_set_ANC(j) = mean(MSPE);        % MSPE for all ensambles for different M
    x_est_mean(j,:) = mean(x_est,1);     % estimated output (ensmable averagr) for different M
end

load 3_3_b_results.mat                   % load pre computed results of MSPE of ALE for various delta aand M 

figure(1); clf
subplot(1,3,1); hold on;
plot(M, MSPE_set(3, :))                  % plot of ALE MSPE for delta=3
plot(M, MSPE_set_ANC)
hold off; box off; grid on;
title('MSPE for ALE and ANC for varying M', 'fontWeight', 'normal')
xlabel('Model order M')
ylabel('MSPE')
legend('ALE, $\Delta$=3', 'ANC', 'location', 'best')
xlim([1,20])

subplot(1,3,2)
plot(x_est_mean(2,:), 'color', [0.8500 0.3250 0.0980])
title('ANC output, $\hat x$, for M=2','Interpreter','latex', 'fontWeight', 'normal')
xlabel('Time [Samples]')
ylabel('ANC output, $\hat x$','Interpreter','latex')
legend('ALE, $\Delta$=3', 'location', 'best')
box off; grid on

% comparing mean square error vs time for ALE and ANC
MSPE_ANC=zeros(num_ens,N);
MSPE_ALE=zeros(num_ens,N);
for i = 1:num_ens
    rng(i)
    wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));    % WGN input.
    eta=filter([1 0 0.5],[1], wgn);                                    % Noise Generation;
    s=x+eta;                                                           % Noisy signal

    [x_est(i,:), w] = ALE_lms(s, 3, u, 4);      % ALE optimal paramters used
    MSPE_ALE(i,:) = ((x_est(i,:) - x).^2);
    
    [x_est(i,:), w] = ANC_lms(s, wgn, u, 2);    % ANC optimal paramters used
    MSPE_ANC(i,:) = ((x_est(i,:) - x).^2);
    
end

subplot(1,3,3)
hold on
plot(10*log10(mean(MSPE_ALE)));
plot(10*log10(mean(MSPE_ANC)));
hold off; box off; grid on
legend('ALE, $\Delta$ = 3, M = 4', 'ANC, M = 2', 'location', 'best')
title('Mean Squared Prediction Error', 'fontWeight', 'normal')
ylabel('MSE, 10*log10( $|x - \hat x |^2$)','Interpreter','latex')
xlabel('Time [Samples]')
ylim([-30, 10])