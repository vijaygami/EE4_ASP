clc;clear;
rng(1)

N=1000;             % number of input samples
p=2;                % AR process order;
u1=0.01;            % adaptation gain
u2=0.05;             
o2=0.25;            % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
b1=[1];             % AR model for sythesis stage
a1=[1, -0.1, -0.8];

%% Evolution of ensamble averaged weights
ave=10;    % Number of ensambles

a1_est_mean=zeros(length(a1)-1,N);
a2_est_mean=zeros(length(a1)-1,N);
for i=1:ave
    rng(i)
    n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input.
    x=filter(b1,a1,n); x=x(501:end, 1);                   % output of filter. Discard first 500 samples due to transient response
    [ y, e1(:,i), a1_est ] = lmsAR( x, u1, p);  
    [ y, e2(:,i), a2_est ] = lmsAR( x, u2, p);   

    a1_est_mean = a1_est_mean + a1_est;
    a2_est_mean = a2_est_mean + a2_est;
end

a1_est_mean=a1_est_mean./ave;
a2_est_mean=a2_est_mean./ave;

figure(1); clf

subplot(1,2,1); hold on
plot([1:N],(a1_est_mean) )
plot([1,N],[-a1(2), -a1(2)])
plot([1,N],[-a1(3), -a1(3)])
xlim([1,N]); ylim([0,0.9]); hold off; box off; grid on
legend('$\hat a_1$','$\hat a_2$', 'True $a_1$', 'True $a_2$', 'location', 'best')
title(['Evolution of Weights, $\mu$ = ', num2str(u1)], 'fontWeight', 'normal')
xlabel('Time [Samples]')
ylabel('Weights')

subplot(1,2,2); hold on
plot([1:N],(a2_est_mean) )
plot([1,N],[-a1(2), -a1(2)])
plot([1,N],[-a1(3), -a1(3)])
xlim([1,N]); ylim([0,0.9]); hold off; box off; grid on
legend('$\hat a_1$','$\hat a_2$', 'True $a_1$', 'True $a_2$', 'location', 'best')
title(['Evolution of Weights, $\mu$ = ', num2str(u2)], 'fontWeight', 'normal')
xlabel('Time [Samples]')
ylabel('Weights')

%% steady state value of coeficients

display('Steady state weights')
ss_u1 = mean(a1_est_mean(:,[800:end]), 2)
ss_u2 = mean(a2_est_mean(:,[300:end]), 2)


%% plots of how coefficient error varies with step size
% WARNING: takes a few minutes to run.(so results are saved and loaded
load '3_1_d_results.mat';
u_set=[0.005:0.0025:0.165]; % set of adaptation gains to find MSE and MSCE for.

% N=5000;
% ave=100;    % Number of ensambles
% 
% Prealocations for speed
% MSCE=zeros(1,N);                % ensamble mean of sum of coeff error squared 
% MSCE_ss=zeros(1,length(u_set)); % steady state MSCE for different adaptation gains
% MSE_ss=zeros(1,length(u_set));  % steady state MSCE for different adaptation gains
% e1=zeros(ave, N);
% 
% for j =1:length(u_set)
%     
%     for i=1:ave
%         rng(i)
%         n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input.
%         x=filter(b1,a1,n); x=x(501:end, 1);                     % output of filter. Discard first 500 samples due to transient response
%         [ y, e1(i,:), a1_est ] = lmsAR( x, u_set(j), p); 
%      
%         MSCE(1,:) = MSCE + (a1_est(1,:)+a1(2)).^2 + (a1_est(2,:)+a1(3)).^2;    % squared coeff error        
%     end
%    
%     MSCE=MSCE./ave;         % Ensamble mean of Squared coef error
%    
%     % Mean squared coefficnet error (in steady state)
%     MSCE_ss(1,j) = mean(MSCE(:,[4000:end]));
%     
%     % Mean squared prediction error
%     error_temp = mean(e1.^2, 1)';
%     MSE_ss(1,j)=10*log10(mean(error_temp(1000:end)));
%     
% end
%%
figure(2); clf
subplot(1,2,1)
plot(u_set, 10*log10(MSCE_ss))
title('Steady State MSCE: 10log10( $|| a - \hat a ||^2$)','Interpreter','latex','fontWeight', 'normal')
xlabel('Adaptation Gain $\mu$','Interpreter','latex')
ylabel('Steady State MSCE','Interpreter','latex','fontWeight', 'normal')
ylim([-30, 10]); xlim([0,0.165]); grid on; box off


subplot(1,2,2)
hold on
plot(u_set, MSE_ss)
plot([0, 0.165],10*log10([0.25,0.25]))
hold off
title('Steady State MSE: 10log10( $|| x - \hat x ||^2$)','Interpreter','latex','fontWeight', 'normal')
xlabel('Adaptation Gain $\mu$','Interpreter','latex')
ylabel('Steady State MSE','Interpreter','latex','fontWeight', 'normal')
leg=legend('Experimental Steady State MSE', 'Minimum Possible MSE', 'location', 'best','Interpreter','latex');  set(leg,'Interpreter','latex')
ylim([-7, -2]); xlim([0,0.165]); grid on; box off
