clc;clear;clf
rng(1)

N=1000;             % number of input samples
p=2;                % AR process order;
u1=0.01;            % adaptation gain
u2=0.05;             
o2=0.25;            % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
b1=[1];             % AR model for sythesis stage
a1=[1, -0.1, -0.8];

gamma=0.1;          % Leaky LMS parameter

%% Evolution of ensamble averaged weights
ave=100;    % Number of ensambles

a1_est_mean=zeros(length(a1)-1,N);
a2_est_mean=zeros(length(a1)-1,N);
for i=1:ave
    rng(i)
    n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input.
    x=filter(b1,a1,n); x=x(501:end, 1);                   % output of filter. Discard first 500 samples due to transient response
    [ y, e1, a1_est ] = leakylmsAR( x, u1, p, gamma);   
    [ y, e2, a2_est ] = leakylmsAR( x, u2, p, gamma);   
    
    a1_est_mean = a1_est_mean + a1_est;
    a2_est_mean = a2_est_mean + a2_est;
end

a1_est_mean=a1_est_mean./ave;
a2_est_mean=a2_est_mean./ave;

figure(1); clf

subplot(1,2,1); hold on
plot([1:N],(a1_est_mean))
plot([1,N],[-a1(2), -a1(2)])
plot([1,N],[-a1(3), -a1(3)])
xlim([1,N]); ylim([0,0.9]); hold off; box off; grid on
title(['Evolution of Weights, $\mu$ = ', num2str(u2), ' $\gamma$ = ', num2str(gamma)], 'fontWeight', 'normal')
xlabel('Iteration')
ylabel('Weights')
legend('$\hat a_1$','$\hat a_2$', 'True $a_1$', 'True $a_2$', 'location', 'best')


subplot(1,2,2); hold on
plot([1:N],(a2_est_mean))
plot([1,N],[-a1(2), -a1(2)])
plot([1,N],[-a1(3), -a1(3)])
xlim([1,N]); ylim([0,0.9]); hold off; box off; grid on
title(['Evolution of Weights, $\mu$ = ', num2str(u2), ' $\gamma$ = ', num2str(gamma)], 'fontWeight', 'normal')
xlabel('Iteration')
ylabel('Weights')
legend('$\hat a_1$','$\hat a_2$', 'True $a_1$', 'True $a_2$', 'location', 'best')


%% steady state value of coeficients

ss_u1 = mean(a1_est_mean(:,[800:end]), 2);
ss_u2 = mean(a2_est_mean(:,[300:end]), 2);



%% plot of steady state value of coeficients as a function of leakage coeficient. slow so load saved results
N=2000;                  % number of input samples
ave=100;                 % Number of ensambles
gamma=[0:0.05:1];        % Leaky LMS parameter
load 3_3_f_results.mat;
% for j=1:length(gamma)
%     
%     a1_est_mean=zeros(length(a1)-1,N);
%     a2_est_mean=zeros(length(a1)-1,N);
%     for i=1:ave
%         rng(i)
%         n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input.
%         x=filter(b1,a1,n); x=x(501:end, 1);                     % output of filter. Discard first 500 samples due to transient response
%         [ y, e1, a1_est ] = leakylmsAR( x, u1, p, gamma(j));   
%         [ y, e2, a2_est ] = leakylmsAR( x, u2, p, gamma(j));   
% 
%         a1_est_mean = a1_est_mean + a1_est;
%         a2_est_mean = a2_est_mean + a2_est;
%     end
% 
%     a1_est_mean=a1_est_mean./ave;
%     a2_est_mean=a2_est_mean./ave;
% 
%     ss(:,j) = [mean(a1_est_mean(:,[1000:end]), 2);mean(a2_est_mean(:,[1000:end]), 2)];
% 
% end

figure(2); clf
subplot(1,2,1); hold on;
plot([0,1],[-a1(2), -a1(2)])
plot([0,1],[-a1(3), -a1(3)])
plot(gamma, ss([1,2],:))
hold off; xlim([0,1]); ylim([0,0.9]); box off; grid on;
legend('True $a_1$', 'True $a_2$','Steady State $\hat a_1$', 'Steady State $\hat a_2$', 'location', 'best')
title('Steady State Weights vs Leakage Coefficient for $\mu$ =0.01', 'fontWeight', 'normal')
xlabel('Leakage Coefficient $\gamma$');
ylabel('Steady State Averaged Weights')

subplot(1,2,2); hold on;

plot([0,1],[-a1(2), -a1(2)])
plot([0,1],[-a1(3), -a1(3)])
plot(gamma, ss([3,4],:))
hold off; xlim([0,1]); ylim([0,0.9]); box off; grid on;
legend('True $a_1$', 'True $a_2$','Steady State $\hat a_1$', 'Steady State $\hat a_2$', 'location', 'best')
title('Steady State Weights vs Leakage Coefficient for $\mu$ =0.05', 'fontWeight', 'normal')
xlabel('Leakage Coefficient $\gamma$');
ylabel('Steady State Averaged Weights')

