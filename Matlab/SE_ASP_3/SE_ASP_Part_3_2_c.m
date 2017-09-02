clc;clear;

N=1000;             % number of input samples
o2=0.5;             % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
q=1;                % MA process order;
b=[1, 0.9];         % MA model for sythesis stage
a=[1];

u_lms=[0.01, 0.05];     % adaptation gain for lms


%% Standard LMS averaged weight errors
ave=100;    % Number of ensambles

b_ben_1=zeros(2,N);
MSCE_ben_1=zeros(1,N);
b_gngd_1=zeros(2,N);
MSCE_gngd_1=zeros(1,N);

for i=1:ave
    
    % MA process synthesis
    rng(i)
    n=randn(N,1); n=n-mean(n); n=n*std./sqrt(var(n));          % WGN input.
    x=filter(b,a,n);                                           % output of filter. Discard first 500 samples due to transient response


    % Benveniste GASS
    [ y, b_est1, e, u] = lmsGASS(n, x, q, 0, 0.025, [], 'Ben');
    b1e=b(1)-b_est1(1,:); b2e=b(2)-b_est1(2,:);     % coeff errors             
    b_ben_1 = b_ben_1 + b_est1;                     % Mean coeficients
    MSCE_ben_1 = MSCE_ben_1 + b1e.^2 + b2e.^2;      % MSCE
    

    % GNGD
    [ y, b_est1, e] = GNGD(n, x, q, 0.1, 0.01, 1);
    b1e=b(1)-b_est1(1,:); b2e=b(2)-b_est1(2,:);     % coeff errors      
    b_gngd_1 = b_gngd_1 + b_est1;                   % Mean coeficients 
    MSCE_gngd_1 = MSCE_gngd_1 + b1e.^2 + b2e.^2;    % MSCE

end

% Ensemble averages
b_ben_1=b_ben_1./ave;
MSCE_ben_1=MSCE_ben_1./ave;
b_gngd_1=b_gngd_1./ave;
MSCE_gngd_1=MSCE_gngd_1./ave;



%% plots
figure(1); clf;
subplot(1,2,1); hold on
plot([1:N], b_ben_1(1,:), 'color',  [0 0.5470 0.7410])
plot([1:N], b_ben_1(2,:), 'color',  [0 0.1470 0.5410])
plot([1:N], b_gngd_1(1,:), 'color',  [0.8500 0.3250 0.0980])
plot([1:N], b_gngd_1(2,:), 'color',  [0.6000 0.2050 0.0580])
hold off; grid on; box off;
title('Evolution of Weights', 'fontWeight', 'normal')
ylabel('Weights')
xlabel('Time [Samples]')
legend('$b_0$ - Benveniste','$b_1$ - Benveniste', '$b_0$ - GNGD', '$b_1$ - GNGD', 'location', 'best')
xlim([0,100]); ylim([0, 1.05])


subplot(1,2,2); hold on
plot([1:N],10*log10(MSCE_ben_1))
plot([1:N],10*log10(MSCE_gngd_1))
hold off; grid on; box off;
title('MSCE: 10log10( $|| b - \hat b ||^2$)','Interpreter','latex','fontWeight', 'normal')
ylabel('10log10( $|| b - \hat b ||$)','Interpreter','latex')
xlabel('Time [Samples]')
legend('Benveniste', 'GNGD', 'location', 'best')
ylim([-350, 10])







