clc;clear;

N=1000;             % number of input samples
o2=0.5;             % Variance of noise
std=sqrt(o2);       % Standard deviation of noise
q=1;                % MA process order;
b=[1, 0.9];         % MA model for sythesis stage
a=[1];

u_lms=[0.01, 0.1];     % adaptation gain for lms


%% Standard LMS averaged weight errors
ave=100;    % Number of ensambles

b_lms_1=zeros(1,N);
b_lms_2=zeros(1,N);
MSCE_lms_1=zeros(1,N);
MSCE_lms_2=zeros(1,N);

b_ben_1=zeros(1,N);
b_ben_2=zeros(1,N);
b_ben_3=zeros(1,N);
b_ben_4=zeros(1,N);
MSCE_ben_1=zeros(1,N);
MSCE_ben_2=zeros(1,N);
MSCE_ben_3=zeros(1,N);
MSCE_ben_4=zeros(1,N);

b_af_1=zeros(1,N);
b_af_2=zeros(1,N);
b_af_3=zeros(1,N);
b_af_4=zeros(1,N);
MSCE_af_1=zeros(1,N);
MSCE_af_2=zeros(1,N);
MSCE_af_3=zeros(1,N);
MSCE_af_4=zeros(1,N);

b_af1_1=zeros(1,N);
b_af1_2=zeros(1,N);
b_af1_3=zeros(1,N);
b_af1_4=zeros(1,N);
MSCE_af1_1=zeros(1,N);
MSCE_af1_2=zeros(1,N);
MSCE_af1_3=zeros(1,N);
MSCE_af1_4=zeros(1,N);

b_mx_1=zeros(1,N);
b_mx_2=zeros(1,N);
b_mx_3=zeros(1,N);
b_mx_4=zeros(1,N);
MSCE_mx_1=zeros(1,N);
MSCE_mx_2=zeros(1,N);
MSCE_mx_3=zeros(1,N);
MSCE_mx_4=zeros(1,N);

for i=1:ave
    
    % MA process synthesis
    rng(i)
    n=randn(N,1); n=n-mean(n); n=n*std./sqrt(var(n));          % WGN input.
    x=filter(b,a,n);                                           % output of filter. Discard first 500 samples due to transient response

    
    % Standard LMS
    [ y, b_est1, e ] = lms(n, x, u_lms(1), q);
    [ y, b_est2, e ] = lms(n, x, u_lms(2), q);
    % Weight errors
    b1e1=b(1)-b_est1(1,:); b2e1=b(2)-b_est1(2,:);              
    b1e2=b(1)-b_est2(1,:); b2e2=b(2)-b_est2(2,:);
    % Mean weight errors of second coeficient only
    b_lms_1 = b_lms_1+b2e1;
    b_lms_2 = b_lms_2+b2e2;
    % MSCE
    MSCE_lms_1 = MSCE_lms_1 + 0.5*(b1e1.^2 + b2e1.^2);
    MSCE_lms_2 = MSCE_lms_2 + 0.5*(b1e2.^2 + b2e2.^2);
    
    
  % Benveniste GASS
    [ y, b_est1, e, u] = lmsGASS(n, x, q, 0, 0.001, [], 'Ben');
    [ y, b_est2, e, u] = lmsGASS(n, x, q, 0, 0.005, [], 'Ben');   
    [ y, b_est3, e, u] = lmsGASS(n, x, q, 0, 0.01, [], 'Ben');
    [ y, b_est4, e, u] = lmsGASS(n, x, q, 0, 0.025, [], 'Ben');
    b1e1=b(1)-b_est1(1,:); b2e1=b(2)-b_est1(2,:);              
    b1e2=b(1)-b_est2(1,:); b2e2=b(2)-b_est2(2,:);
    b1e3=b(1)-b_est3(1,:); b2e3=b(2)-b_est3(2,:);
    b1e4=b(1)-b_est4(1,:); b2e4=b(2)-b_est4(2,:);
    % Mean weight errors of second coeficient only
    b_ben_1 = b_ben_1+b2e1;
    b_ben_2 = b_ben_2+b2e2;
    b_ben_3 = b_ben_4+b2e3;
    b_ben_4 = b_ben_4+b2e4;
    % MSCE
    MSCE_ben_1 = MSCE_ben_1 + 0.5*(b1e1.^2 + b2e1.^2);
    MSCE_ben_2 = MSCE_ben_2 + 0.5*(b1e2.^2 + b2e2.^2);
    MSCE_ben_3 = MSCE_ben_3 + 0.5*(b1e3.^2 + b2e3.^2);
    MSCE_ben_4 = MSCE_ben_4 + 0.5*(b1e4.^2 + b2e4.^2);   

    
  % Ang & Farhang GASS varying rho only
    [ y, b_est1, e, u] = lmsGASS(n, x, q, 0, 0.001, 0.5, 'A&F');  
    [ y, b_est2, e, u] = lmsGASS(n, x, q, 0, 0.01, 0.5, 'A&F');    
    [ y, b_est3, e, u] = lmsGASS(n, x, q, 0, 0.025, 0.5, 'A&F');  
    [ y, b_est4, e, u] = lmsGASS(n, x, q, 0, 0.05, 0.5, 'A&F'); 
    b1e1=b(1)-b_est1(1,:); b2e1=b(2)-b_est1(2,:);              
    b1e2=b(1)-b_est2(1,:); b2e2=b(2)-b_est2(2,:);
    b1e3=b(1)-b_est3(1,:); b2e3=b(2)-b_est3(2,:);
    b1e4=b(1)-b_est4(1,:); b2e4=b(2)-b_est4(2,:);
    % Mean weight errors of second coeficient only
    b_af_1 = b_af_1+b2e1;
    b_af_2 = b_af_2+b2e2;
    b_af_3 = b_af_3+b2e3;
    b_af_4 = b_af_4+b2e4;
    % MSCE
    MSCE_af_1 = MSCE_af_1 + 0.5*(b1e1.^2 + b2e1.^2);
    MSCE_af_2 = MSCE_af_2 + 0.5*(b1e2.^2 + b2e2.^2);   
    MSCE_af_3 = MSCE_af_3 + 0.5*(b1e3.^2 + b2e3.^2);
    MSCE_af_4 = MSCE_af_4 + 0.5*(b1e4.^2 + b2e4.^2);
    
    
  % Ang & Farhang GASS varying alpha only
    [ y, b_est1, e, u] = lmsGASS(n, x, q, 0, 0.025, 0.1, 'A&F');  
    [ y, b_est2, e, u] = lmsGASS(n, x, q, 0, 0.025, 0.4, 'A&F');    
    [ y, b_est3, e, u] = lmsGASS(n, x, q, 0, 0.025, 0.7, 'A&F');  
    [ y, b_est4, e, u] = lmsGASS(n, x, q, 0, 0.025, 0.95, 'A&F'); 

    b1e1=b(1)-b_est1(1,:); b2e1=b(2)-b_est1(2,:);              
    b1e2=b(1)-b_est2(1,:); b2e2=b(2)-b_est2(2,:);
    b1e3=b(1)-b_est3(1,:); b2e3=b(2)-b_est3(2,:);
    b1e4=b(1)-b_est4(1,:); b2e4=b(2)-b_est4(2,:);
    % Mean weight errors of second coeficient only
    b_af1_1 = b_af1_1+b2e1;
    b_af1_2 = b_af1_2+b2e2;
    b_af1_3 = b_af1_3+b2e3;
    b_af1_4 = b_af1_4+b2e4;
    % MSCE
    MSCE_af1_1 = MSCE_af1_1 + 0.5*(b1e1.^2 + b2e1.^2);
    MSCE_af1_2 = MSCE_af1_2 + 0.5*(b1e2.^2 + b2e2.^2);   
    MSCE_af1_3 = MSCE_af1_3 + 0.5*(b1e3.^2 + b2e3.^2);
    MSCE_af1_4 = MSCE_af1_4 + 0.5*(b1e4.^2 + b2e4.^2);

    
    % Matthews and Xie GASS
    [ y, b_est1, e, u] = lmsGASS(n, x, q, 0, 0.001, [], 'M&X');  
    [ y, b_est2, e, u] = lmsGASS(n, x, q, 0, 0.002, [], 'M&X');  
    [ y, b_est3, e, u] = lmsGASS(n, x, q, 0, 0.005, [], 'M&X');  
    [ y, b_est4, e, u] = lmsGASS(n, x, q, 0, 0.01, [], 'M&X'); 
    b1e1=b(1)-b_est1(1,:); b2e1=b(2)-b_est1(2,:);              
    b1e2=b(1)-b_est2(1,:); b2e2=b(2)-b_est2(2,:);
    b1e3=b(1)-b_est3(1,:); b2e3=b(2)-b_est3(2,:);              
    b1e4=b(1)-b_est4(1,:); b2e4=b(2)-b_est4(2,:);
    % Mean weight errors of second coeficient only
    b_mx_1 = b_mx_1+b2e1;
    b_mx_2 = b_mx_2+b2e2;
    b_mx_3 = b_mx_3+b2e3;
    b_mx_4 = b_mx_4+b2e4;
    % MSCE
    MSCE_mx_1 = MSCE_mx_1 + 0.5*(b1e1.^2 + b2e1.^2);
    MSCE_mx_2 = MSCE_mx_2 + 0.5*(b1e2.^2 + b2e2.^2);
    MSCE_mx_3 = MSCE_mx_3 + 0.5*(b1e3.^2 + b2e3.^2);
    MSCE_mx_4 = MSCE_mx_4 + 0.5*(b1e4.^2 + b2e4.^2);   
end

% Ensemble averages
b_lms_1=b_lms_1./ave;
b_lms_2=b_lms_2./ave;
MSCE_lms_1=MSCE_lms_1./ave;
MSCE_lms_2=MSCE_lms_2./ave;

b_ben_1=b_ben_1./ave;
b_ben_2=b_ben_2./ave;
b_ben_3=b_ben_3./ave;
b_ben_4=b_ben_4./ave;
MSCE_ben_1=MSCE_ben_1./ave;
MSCE_ben_2=MSCE_ben_2./ave;
MSCE_ben_3=MSCE_ben_3./ave;
MSCE_ben_4=MSCE_ben_4./ave;

b_af_1=b_af_1./ave;
b_af_2=b_af_2./ave;
b_af_3=b_af_3./ave;
b_af_4=b_af_4./ave;
MSCE_af_1=MSCE_af_1./ave;
MSCE_af_2=MSCE_af_2./ave;
MSCE_af_3=MSCE_af_3./ave;
MSCE_af_4=MSCE_af_4./ave;

b_af1_1=b_af1_1./ave;
b_af1_2=b_af1_2./ave;
b_af1_3=b_af1_3./ave;
b_af1_4=b_af1_4./ave;
MSCE_af1_1=MSCE_af1_1./ave;
MSCE_af1_2=MSCE_af1_2./ave;
MSCE_af1_3=MSCE_af1_3./ave;
MSCE_af1_4=MSCE_af1_4./ave;

b_mx_1=b_mx_1./ave;
b_mx_2=b_mx_2./ave;
b_mx_3=b_mx_3./ave;
b_mx_4=b_mx_4./ave;
MSCE_mx_1=MSCE_mx_1./ave;
MSCE_mx_2=MSCE_mx_2./ave;
MSCE_mx_3=MSCE_mx_3./ave;
MSCE_mx_4=MSCE_mx_4./ave;
%% plots of each GASS algo wth varying parameters

figure(1); clf;
subplot(4,2,1); hold on
plot(b_ben_1)
plot(b_ben_2)
plot(b_ben_3)
plot(b_ben_4)
xlim([1,100]); hold off; grid on
legend('$\rho$=0.001', '$\rho$=0.005','$\rho$=0.01', '$\rho$=0.025', 'location', 'best')
title('Benveniste Weight Error', 'fontWeight', 'normal')
ylabel('Weight Error, $b1$ - $\hat b1$')
subplot(4,2,2); hold on
plot([1:N],10*log10(MSCE_ben_1))
plot([1:N],10*log10(MSCE_ben_2))
plot([1:N],10*log10(MSCE_ben_3))
plot([1:N],10*log10(MSCE_ben_4))
xlim([1,1000]); hold off; grid on; box off
legend('$\rho$=0.001', '$\rho$=0.005','$\rho$=0.01', '$\rho$=0.025', 'location', 'best')
title('Benveniste MSCE: 10log10( $|| b - \hat b ||^2$)','fontWeight', 'normal')
ylabel('10log10( $||b - \hat b||^2$)')
ylim([-400, 75])


subplot(4,2,3); hold on
plot(b_af_1)
plot(b_af_2)
plot(b_af_3)
plot(b_af_4)
xlim([1,100]); hold off; grid on
legend('$\rho$=0.001', '$\rho$=0.01','$\rho$=0.025', '$\rho$=0.05', 'location', 'best')
title('Ang $\&$ Farhang (fixed $\alpha$ = 0.5) Weight  Error ', 'fontWeight', 'normal')
ylabel('Weight Error, $b1$ - $\hat b1$')
subplot(4,2,4); hold on
plot([1:N],10*log10(MSCE_af_1))
plot([1:N],10*log10(MSCE_af_2))
plot([1:N],10*log10(MSCE_af_3))
plot([1:N],10*log10(MSCE_af_4))
xlim([1,1000]); hold off; grid on; box off
legend('$\rho$=0.001', '$\rho$=0.01','$\rho$=0.025', '$\rho$=0.05', 'location', 'best')
title('Ang $\&$ Farhang (fixed $\alpha$ = 0.5) MSCE : 10log10( $|| b - \hat b ||^2$)','Interpreter','latex','fontWeight', 'normal')
ylabel('10log10( $||b - \hat b||^2$)')
ylim([-400, 75])


subplot(4,2,5); hold on
plot(b_af1_1)
plot(b_af1_2)
plot(b_af1_3)
plot(b_af1_4)
xlim([1,100]); hold off; grid on
legend('$\alpha$=0.1', '$\alpha$=0.4','$\alpha$=0.7', '$\alpha$=0.95', 'location', 'best')
title('Ang $\&$ Farhang (fixed $\rho$ = 0.025) Weight Error', 'fontWeight', 'normal')
ylabel('Weight Error, $b1$ - $\hat b1$')
subplot(4,2,6); hold on
plot([1:N],10*log10(MSCE_af1_1))
plot([1:N],10*log10(MSCE_af1_2))
plot([1:N],10*log10(MSCE_af1_3))
plot([1:N],10*log10(MSCE_af1_4))
xlim([1,1000]); hold off; grid on; box off
legend('$\alpha$=0.1', '$\alpha$=0.4','$\alpha$=0.7', '$\alpha$=0.95', 'location', 'best')
title('Ang \& Farhang (fixed $\rho$ = 0.025) MSCE: 10log10( $|| b - \hat b ||^2$)','Interpreter','latex','fontWeight', 'normal')
ylabel('10log10( $||b - \hat b||^2$)','Interpreter','latex')
ylim([-400, 75])


subplot(4,2,7); hold on
plot(b_mx_1)
plot(b_mx_2)
plot(b_mx_3)
plot(b_mx_4)
xlim([1,100]); hold off; grid on
legend('$\rho$=0.001', '$\rho$=0.002','$\rho$=0.005', '$\rho$= 0.01', 'location', 'best')
title('Matthews $\&$ Xie Weight Error', 'fontWeight', 'normal')
ylabel('Weight Error, $b1$ - $\hat b1$', 'Interpreter','latex')
subplot(4,2,8); hold on
plot([1:N],10*log10(MSCE_mx_1))
plot([1:N],10*log10(MSCE_mx_2))
plot([1:N],10*log10(MSCE_mx_3))
plot([1:N],10*log10(MSCE_mx_4))
xlim([1,1000]); hold off; grid on; box off
legend('$\rho$=0.001', '$\rho$=0.002','$\rho$=0.005', '$\rho$= 0.01', 'location', 'best')
title('Matthews \& Xie MSCE: 10log10( $|| b - \hat b ||^2$)', 'Interpreter','latex','fontWeight', 'normal')
ylabel('10log10( $||b - \hat b||^2$)','Interpreter','latex')
ylim([-400, 75])

for i=1:8
    subplot(4,2,i)
    xlabel('Time [Samples]')
end

%% plots of three GASS algos with optimal parameters compared to standard LMS
figure(2); clf
subplot(1,2,1); hold on
plot(b_ben_4)
plot(b_af1_4)
plot(b_mx_4)
plot(b_lms_1)
plot(b_lms_2)
xlim([1,150]); hold off; grid on
legend('Benveniste $\rho$=0.025', 'Ang $\&$ Farhang $\rho$=0.025 $\alpha$=0.95','Matthews $\&$ Xie $\rho$=0.01', 'LMS $\mu$=0.01 ', 'LMS $\mu$=0.1 ', 'location', 'best')
title('GASS Algorithms vs LMS, Weight Error', 'fontWeight', 'normal')
ylabel('Weight Error, $b1$ - $\hat b1$', 'Interpreter','latex')
xlabel('Time [Samples]')
ylim([0,1.05])
subplot(1,2,2); hold on
plot([1:N],10*log10(MSCE_ben_4))
plot([1:N],10*log10(MSCE_af1_4))
plot([1:N],10*log10(MSCE_mx_4))
plot([1:N],10*log10(MSCE_lms_1))
plot([1:N],10*log10(MSCE_lms_2))
xlim([1,1000]); hold off; grid on; box off
legend('Benveniste $\rho$=0.025', 'Ang $\&$ Farhang $\rho$=0.025 $\alpha$=0.95','Matthews $\&$ Xie $\rho$=0.01', 'LMS $\mu$=0.01 ', 'LMS $\mu$=0.1 ', 'location', 'best')
title('GASS Algorithms vs LMS, MSCE: 10log10( $|| b - \hat b ||^2$)','Interpreter','latex','fontWeight', 'normal')
ylabel('10log10( $|| b - \hat b ||^2$)','Interpreter','latex')
xlabel('Time [Samples]')
ylim([-350, 50])