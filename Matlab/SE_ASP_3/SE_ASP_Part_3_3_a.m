clc;clear;

N=1000;             % Number of input samples
wo=0.01*pi;         % Frequency of input signal
o2=1;               % Variance of white noise
std=sqrt(o2);       % Standard deviation of the white noise

M=5;                % ALE line enchancer filter  order;
delta=[1, 2, 3, 6];

u=0.01;             % Filter adaptation gain
num_ens=100;        % number of independant trials

% Signal generation
x=sin(wo.*[1:N]);


%% Effects of Delta on the ALE

for k=1:length(delta)
    for i = 1:num_ens
        rng(i)
        wgn=randn(1,N); wgn=wgn-mean(wgn); wgn=wgn*std./sqrt(var(wgn));  % WGN input.
        eta=filter([1 0 0.5],[1], wgn);                                  % Noise Generation
        s=x+eta;                                                         % Noisy signal    
        [x_est(i,:), e, w] = ALE_lms(s, delta(k), u, M); % LMS algorithm
    end
        x_est_m(k,:)=mean(x_est);
        e_m(k,:)=mean((x_est-ones(num_ens, 1)*x).^2);
        MSPE(k,1)=mean(e_m(k,:));

end

%% Plots 

figure(1); clf;
for i=1:length(delta)
    subplot(2,length(delta), i);
    plot([1:N], x_est_m(i,:))
    xlabel('Time [Samples]')
    ylabel('Predicted Output, $\hat x$','Interpreter','latex')
    title(['Prediction, $\Delta$=',num2str(delta(i)), ', M=', num2str(M)], 'fontWeight', 'normal')
    ylim([-1.1,1.2]); grid on; box off;
end

for i=1:length(delta)
    subplot(2,length(delta), i+length(delta)); hold on
    plot([1:N], e_m(i,:))
    text(200,0.9,['MSPE = ', num2str(MSPE(i))])
    hold off
    xlabel('Time [Samples]')
    ylabel('MSE, $(x - \hat x$)$^2$','Interpreter','latex')
    title(['Prediction Error, $\Delta$=',num2str(delta(i)), ', M=', num2str(M)], 'fontWeight', 'normal')
    grid on; box off;
    ylim([0,1])
end

