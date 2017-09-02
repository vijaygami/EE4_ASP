clc;clear;

N=1000;             % Number of input samples
wo=0.01*pi;         % Frequency of input signal
o2=1;               % Variance of white noise
std=sqrt(o2);       % Standard deviation of the white noise

M=[1:1:25];         % ALE line enchancer filter  order;
delta=[1:1:120];

u=0.01;             % Filter adaptation gain
num_ens=100;        % number of independant trials

% Signal generation
x=sin(wo.*[1:N]);

%% Effects on MSPE by varying M and Delta 
%the results are saved as a file since this computation takes very long (10 min)

% x_est=zeros(num_ens,N); % pre-allocate memory for speed
% MSPE=zeros(num_ens,1);
% MSPE_set=zeros(length(delta),length(M));
% for k=1:length(delta)
%     for j=1:length(M)
%         for i = 1:num_ens
%             rng(i)
%             % Noise Generation;
%             wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));  % WGN input.
%             eta=filter([1 0 0.5],[1], wgn);
%             s=x+eta; % Noisy signal
%             
%             [x_est(i,:), e, w] = ALE_lms(s, delta(k), u, M(j)); % LMS algorithm
%             MSPE(i,1) = mean((x_est(i,:) - x).^2);              % Mean square prediciton error
%         end
%             MSPE_set(k,j)=mean(MSPE);                             % Mean MSPE of all ensambles
%     end
%     k %progress counter since it is slow
% end
% save('3_3_b_results', 'MSPE_set')       % Save variable for future use since computation is long

% plots
load 3_3_b_results.mat

figure(1); clf
subplot(1,2,1); hold on
plot(delta, MSPE_set(:, 5))
plot(delta, MSPE_set(:, 10))
plot(delta, MSPE_set(:, 15))
plot(delta, MSPE_set(:, 20))
hold off; box off; grid on
title('MSPE as a function of $\Delta$','fontWeight', 'normal')
xlabel('Delay, $\Delta$')
ylabel('MSPE')
legend('M = 5', 'M = 10', 'M = 15', 'M = 10', 'location', 'best' )

subplot(1,2,2); hold on
plot(M, MSPE_set(2, :))
plot(M, MSPE_set(3, :))
plot(M, MSPE_set(10, :))
plot(M, MSPE_set(50, :))
plot(M, MSPE_set(100, :))

hold off; box off; grid on
title('MSPE as a function of M','fontWeight', 'normal')
xlabel('Model order, M')
ylabel('MSPE')
ylim([0.25, 0.85])
legend('$\Delta$ = 2', '$\Delta$ = 3', '$\Delta$ = 10', '$\Delta$ = 50', '$\Delta$ = 100', 'location', 'best' )


%% plots of prediction when varying M

M=[5 10 15 20];         % ALE line enchancer filter  order;
delta=3;

for k=1:length(M)
    for i = 1:num_ens
        rng(i)
        wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));  % WGN input.
        eta=filter([1 0 0.5],[1], wgn);                                  % Noise Generation
        s=x+eta;                                                         % Noisy signal    
        [x_est(i,:), e, w] = ALE_lms(s, delta, u, M(k)); % LMS algorithm
    end
        x_est_m(k,:)=mean(x_est);
        e_m(k,:)=mean((x_est-ones(num_ens, 1)*x).^2);
        MSPE(k,1)=mean(e_m(k,:));

end


% plots
figure(2); clf;
for i=1:length(M)
    subplot(2,length(M), i);
    plot([1:N], x_est_m(i,:))
    xlabel('Time [Samples]')
    ylabel('Predicted Output, $\hat x$','Interpreter','latex')
    title(['Prediction, $\Delta$=',num2str(delta), ', M=', num2str(M(i))], 'fontWeight', 'normal')
    grid on; box off;
    ylim([-1.1,1.2])

end

for i=1:length(M)
    subplot(2,length(M), i+length(M));
    plot([1:N], e_m(i,:)); hold on
    text(200,1.05,['MSPE = ', num2str(MSPE(i))])
    hold off
    xlabel('Time [Samples]')
    ylabel('MSE, $(x - \hat x$)$^2$','Interpreter','latex')
    title(['Prediction Error, $\Delta$=',num2str(delta), ', M=', num2str(M(i))], 'fontWeight', 'normal')
    grid on; box off;
    ylim([0,1.2])
end
