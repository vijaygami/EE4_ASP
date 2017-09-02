clc;clear
rng(1)

N=1000;             % number of input samples
p=2;                % AR process order;
u=[0.01 0.05];      % adaptation gain
o2=0.25;            % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
b1=[1];             % AR model for sythesis stage
a1=[1, -0.1, -0.8];

%% synthesis stage    
n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input, zero mean and variance = o2
x=filter(b1,a1,n); x=x([501:end], 1);                   % Output of filter. Discard first 500 samples due to transient response

%% LMS algorithm and plot of error power

for i = 1:length(u)
    [ y, e(i,:), a ] = lmsAR( x, u(i), p);
    legendInfo{i} = ['$\mu = $' num2str(u(i))];             % legend entries for plot
end
figure(1); clf
subplot(2,1,1); hold on
plot([1:N], 10*log10(e.^2));
title('Learning Curve for AR(2), x(n) = 0.1x(n-1) + 0.8x(n-2) + $\eta $(n)', 'FontWeight', 'normal');
line([1, N],[10*log10(o2), 10*log10(o2)], 'lineWidth', 1.5, 'color', 'black') % theoretical Min MSE

xlabel('Time [Samples]')
ylabel('Squared Predicition Error [dB]');
legend ([legendInfo, 'Theoretical Minimum MSE'])
ylim([-60, 25]); hold off; box off; grid on

%% plot of mean square error

n_ave=100;    % Number of averages to perform
MSE=zeros(2,N);  % Will store MSE

for j=1:n_ave  
    % synthesis stage
    rng(j)   
    n=randn(N+500,1); n=n-mean(n); n=n*std./sqrt(var(n));   % WGN input, zero mean and variance = o2
    x=filter(b1,a1,n); x=x(501:end, 1);                     % output of filter. Discard first 500 samples due to transient response

    for i = 1:length(u)
        [ y, e(i,:), a ] = lmsAR( x, u(i), p);   
        e(i,:)=e(i,:).^2;
        MSE(i,:)=MSE(i,:)+e(i,:);
    end
    
end

MSE=MSE./n_ave;  % Mean of error squared

% Plots
subplot(2,1,2)
hold on
plot([1:N], 10*log10(MSE))
line([1, N],[10*log10(o2), 10*log10(o2)], 'lineWidth', 1.5, 'color', 'black') % theoretical Min MSE
legend('$\mu $= 0.01','$\mu $= 0.05', 'Theoretical Minimum MSE')
hold off
xlabel('Time [Samples]');
ylabel('MSE [dB]');
title('Averaged Learning Curve for AR(2), x(n) = 0.1x(n-1) + 0.8x(n-2) + $\eta $(n). 100 Averages', 'FontWeight', 'normal');
ylim([-8, 1]); box off; grid on




