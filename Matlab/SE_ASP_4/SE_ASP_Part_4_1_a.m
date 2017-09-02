clc;clear;

N=2000;                     % number of input samples
q=1;                        % MA process order
u=0.025;                    % adaptation gain         
o2=1;                       % Variance of noise
std=sqrt(o2);               % Standard deviation of noise
b=[1.5+1i, 2.5-0.5i];       % array of MA coef for sythesis stage

%% Evolution of ensamble averaged weights
num_ens=100;        % Number of ensambles

% a1_est_mean=zeros(length(a)-1,N+1);
e_aclms = zeros(num_ens, N);
e_clms = zeros(num_ens, N);
y_est_aclms = zeros(1, N);
y_est_clms = zeros(1, N);

h_est_clms_ave = zeros(q+1,N);
h_est_aclms_ave = zeros(q+1,N);
g_est_aclms_ave = zeros(q+1,N);

for i=1:num_ens
    rng(i);
    cwgn = wgn(1,N,0,'complex');
    y=zeros(1,N);
    
    for k=2:N
        y(k) = cwgn(k) + b(1)*cwgn(k-1) + b(2)*conj(cwgn(k-1));
    end 

    [y_est_clms(i,:), h_est_clms, e_clms(i,:)] = CLMS(cwgn, y, u, q);  
    [y_est_aclms(i,:), h_est_aclms, g_est_aclms, e_aclms(i,:)] = ACLMS(cwgn, y, u, q);  

     h_est_clms_ave = h_est_clms_ave + h_est_clms;
     h_est_aclms_ave = h_est_aclms_ave + h_est_aclms;
     g_est_aclms_ave = g_est_aclms_ave + g_est_aclms;    
end


h_est_clms_ave = h_est_clms_ave./num_ens;
h_est_aclms_ave = h_est_aclms_ave./num_ens;
g_est_aclms_ave = g_est_aclms_ave./num_ens;

e_clms=mean(abs(e_clms).^2, 1);
e_aclms=mean(abs(e_aclms).^2, 1);

%% plots of circularity of the data

figure(1); clf; hold on
plot(real(y), imag(y), '.')
plot([-12, 12], [0, 0], '--', 'color', 'red')
plot([0, 0], [-12, 12], '--', 'color', 'red')
hold off
c=abs(circ_coef(y));
axis([-12, 12, -12, 12])
title(['Circularity Diagram, Circ Coeff = ', num2str(c)], 'fontWeight', 'normal')
xlabel('Real Part')
ylabel('Imaginary part')


%% plots to compare clms and aclms

figure(2); clf; % CLMS
subplot(2,2,1)
plot(10*log10(e_clms)); hold on;
plot([0, 2000], [10*log10(mean(e_clms(1000:end))), 10*log10(mean(e_clms(1000:end)))], '--', 'color', 'red', 'lineWidth', 1)
hold off; box on; grid on;
title(['CLMS Learning curve of a WLMA(1) Process'], 'fontWeight', 'normal', 'Interpreter','latex')
xlabel('Time [Samples]', 'Interpreter','latex')
ylabel('MSE, 10log10( $|y - \hat y|^2$)','Interpreter','latex','fontWeight', 'normal')
legend('MSE', 'Steady State MSE')
ylim([0,15])

subplot(2,2,2); hold on
plot(real(h_est_clms_ave).', 'lineWidth', 1)
plot(imag(h_est_clms_ave).', 'lineWidth', 1)
hold off; box on; grid on;
legend('$\Re(h_1)$', '$\Re(h_2)$', '$\Im(h_1)$', '$\Im(h_2)$')
xlim([0,400])
title(['Evolution of Weights for CLMS'], 'fontWeight', 'normal', 'Interpreter','latex')
xlabel('Time [Samples]', 'Interpreter','latex')
ylabel('Weights, h','Interpreter','latex','fontWeight', 'normal', 'Interpreter','latex')


%ACLMS
subplot(2,2,3);
plot(10*log10(e_aclms)); hold on;
plot([0, 2000], [10*log10(mean(e_aclms(1800:end))), 10*log10(mean(e_aclms(1800:end)))], '--', 'color', 'red', 'lineWidth', 1)
hold off; box on; grid on;
title(['ACLMS Learning curve of a WLMA(1) Process'], 'fontWeight', 'normal', 'Interpreter','latex')
xlabel('Time [Samples]', 'Interpreter','latex')
ylabel('MSE, 10log10( $|y - \hat y|^2$)','Interpreter','latex','fontWeight', 'normal')
legend('MSE', 'Steady State MSE')

subplot(2,2,4); hold on
plot(real(h_est_aclms_ave).', 'lineWidth', 1)
plot(imag(h_est_aclms_ave).', 'lineWidth', 1)
plot(real(g_est_aclms_ave).', 'lineWidth', 1)
plot(imag(g_est_aclms_ave).', 'lineWidth', 1)
hold off; box on; grid on;
xlim([0,400]); ylim([-1.5, 3])
legend('$\Re(h_1)$', '$\Re(h_2)$', '$\Im(h_1)$', '$\Im(h_2)$', '$\Re(g_1)$', '$\Re(g_2)$', '$\Im(g_1)$', '$\Im(g_2)$')
title(['Evolution of Weights for ACLMS'], 'fontWeight', 'normal', 'Interpreter','latex')
xlabel('Time [Samples]', 'Interpreter','latex')
ylabel('Weights, h, g','Interpreter','latex','fontWeight', 'normal', 'Interpreter','latex')