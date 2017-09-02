clc;clear;
% load data as complex vectors
load('high-wind.mat'); v_h = v_east + 1i*v_north;
load('medium-wind.mat'); v_m = v_east + 1i*v_north;
load('low-wind.mat'); v_l = v_east + 1i*v_north;
N=length(v_h);

% preprocessing:
v_h=detrend(v_h);
v_m=detrend(v_m);
v_l=detrend(v_l);

%% Circularity plots of wind speed data

figure(1); clf;
subplot(1,3,1); hold on
plot(real(v_l), imag(v_l), '.', 'markerSize', 6)
rectangle('Position',[-1 -1 2 2]*0.001 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 2)
rectangle('Position',[-1 -1 2 2]*0.1 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*0.2 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*0.3 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*0.4 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*0.5 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*0.6 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
hold off; xlim([-0.45, 0.45]), ylim([-0.45, 0.45])
title({'Circularity plot, Low Wind Speed',['Circularity Coefficient = ', num2str(abs(circ_coef(v_l))), ' $\angle$', num2str(angle(circ_coef(v_l)))]}, 'fontWeight', 'normal')


subplot(1,3,2);
plot(real(v_m), imag(v_m), '.', 'markerSize', 6); hold on
rectangle('Position',[-1 -1 2 2]*0.005 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 2)
rectangle('Position',[-1 -1 2 2]*0.5 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*1 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*1.5 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*2 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*2.5 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*3 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
hold off; xlim([-2, 2]), ylim([-2, 2])
title({'Circularity plot, Medium Wind Speed',['Circularity Coefficient = ', num2str(abs(circ_coef(v_m))), ' $\angle$', num2str(angle(circ_coef(v_m)))]}, 'fontWeight', 'normal')


subplot(1,3,3); hold on
plot(real(v_h), imag(v_h), '.', 'markerSize', 6)
rectangle('Position',[-1 -1 2 2]*0.005 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 2)
rectangle('Position',[-1 -1 2 2]*1 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*2 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*3 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*4 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*5 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
rectangle('Position',[-1 -1 2 2]*6 ,'Curvature',[1 1], 'lineStyle', '--', 'edgeColor', 'red', 'lineWidth', 0.5)
hold off; xlim([-4.5, 4.5]), ylim([-4.5, 4.5])
title({'Circularity plot, High Wind Speed',['Circularity Coefficient = ', num2str(abs(circ_coef(v_h))), ' $\angle$', num2str(angle(circ_coef(v_h)))]}, 'fontWeight', 'normal')

% common parts of all figures in a loop
for i=1:3
    subplot(1,3,i)  
    xlabel('Real Part')
    ylabel('Imaginary part')
    hold on
   % plot([-10, 10], [0, 0], '--', 'color', 'red')
   % plot([0, 0], [-10, 10], '--', 'color', 'red')
    hold off; box off, grid off;
end

%% Prediction errors for ACLMS and CLMS
p=[1:1:30];

ul=0.002;    % adaptation gains for different wind speeds
um=0.001;
uh=0.00025;

for i = 1: length(p)
  % CLMS
    [y_est_clms_l, ~, ~ ] = CLMS_AR(v_l, ul, p(i));
    mspe_clms_low(1,i) = mean(abs(v_l-y_est_clms_l.').^2);

    [y_est_clms_m, ~, ~ ] = CLMS_AR(v_m, um, p(i));
    mspe_clms_medium(1,i) = mean(abs(v_m-y_est_clms_m.').^2);
    
    [y_est_clms_h, ~, ~ ] = CLMS_AR(v_h, uh, p(i));
    mspe_clms_high(1,i) = mean(abs(v_h-y_est_clms_h.').^2);
    
    
  % ACLMS
    [y_est_aclms_l, ~, ~, ~ ] = ACLMS_AR(v_l, ul, p(i));
    mspe_aclms_low(1,i) = mean(abs(v_l-y_est_aclms_l.').^2);

    [y_est_aclms_m, ~, ~, ~ ] = ACLMS_AR(v_m, um, p(i));
    mspe_aclms_medium(1,i) = mean(abs(v_m-y_est_aclms_m.').^2);
    
    [y_est_aclms_h, ~, ~, ~ ] = ACLMS_AR(v_h, uh, p(i));
    mspe_aclms_high(1,i) = mean(abs(v_h-y_est_aclms_h.').^2);
    
end

%% Plots of prediction errors for ACLMS and CLMS
figure(2); clf
subplot(1,3,1); hold on
plot(p, 10*log10(mspe_clms_low))
plot(p, 10*log10(mspe_aclms_low))
hold off; %ylim([-20.5, -17.5])
title('MSPE for Low Wind Speed', 'fontWeight', 'normal')

subplot(1,3,2); hold on
plot(p, 10*log10(mspe_clms_medium))
plot(p, 10*log10(mspe_aclms_medium))
hold off; %ylim([-11, -8])
title('MSPE for Medium Wind Speed', 'fontWeight', 'normal')

subplot(1,3,3); hold on
plot(p, 10*log10(mspe_clms_high))
plot(p, 10*log10(mspe_aclms_high))
hold off; %ylim([-5, -2])
title('MSPE for High Wind Speed', 'fontWeight', 'normal')

for i=1:3
    subplot(1,3,i)
    legend('CLMS','ACLMS')
    xlabel('Model Order')
    ylabel('MSPE')
    box off; grid on
end

%% Circularity plots showing actual data and prediction
% 
% p=[10];
% for i = 1: length(p)
%   % CLMS
%     [y_est_clms_l, ~, ~ ] = CLMS_AR(v_l, ul, p(i));
%     mspe_clms_low(1,i) = mean(abs(v_l-y_est_clms_l.').^2);
% 
%     [y_est_clms_m, ~, ~ ] = CLMS_AR(v_m, um, p(i));
%     mspe_clms_medium(1,i) = mean(abs(v_m-y_est_clms_m.').^2);
%     
%     [y_est_clms_h, ~, ~ ] = CLMS_AR(v_h, uh, p(i));
%     mspe_clms_high(1,i) = mean(abs(v_h-y_est_clms_h.').^2);
%     
%     
%   % ACLMS
%     [y_est_aclms_l, ~, ~, ~ ] = ACLMS_AR(v_l, ul, p(i));
%     mspe_aclms_low(1,i) = mean(abs(v_l-y_est_aclms_l.').^2);
% 
%     [y_est_aclms_m, ~, ~, ~ ] = ACLMS_AR(v_m, um, p(i));
%     mspe_aclms_medium(1,i) = mean(abs(v_m-y_est_aclms_m.').^2);
%     
%     [y_est_aclms_h, ~, ~, ~ ] = ACLMS_AR(v_h, uh, p(i));
%     mspe_aclms_high(1,i) = mean(abs(v_h-y_est_aclms_h.').^2);
%     
% end
% 
% figure(3); clf;
% subplot(1,3,1); hold on;
% plot(real(v_l), imag(v_l), '.', 'markerSize', 3)
% plot(real(y_est_clms_l), imag(y_est_clms_l), '.', 'markerSize', 3)
% plot(real(y_est_aclms_l), imag(y_est_aclms_l), '.', 'markerSize', 3)
% plot([-10, 10], [0, 0], '--', 'color', 'red')
% plot([0, 0], [-10, 10], '--', 'color', 'red')
% hold off; xlim([-0.45, 0.45]), ylim([-0.45, 0.45])
% title('Low Wind Speed', 'fontWeight', 'normal')
% legend('Data', 'CLMS', 'ACLMS')
% 
% subplot(1,3,2); hold on
% plot(real(v_m), imag(v_m), '.', 'markerSize', 3)
% plot(real(y_est_clms_h), imag(y_est_clms_h), '.', 'markerSize', 3)
% plot(real(y_est_aclms_h), imag(y_est_aclms_h), '.', 'markerSize', 3)
% plot([-10, 10], [0, 0], '--', 'color', 'red')
% plot([0, 0], [-10, 10], '--', 'color', 'red')
% hold off; xlim([-2, 2]), ylim([-2, 2])
% title('Medium Wind Speed', 'fontWeight', 'normal')
% legend('Data', 'CLMS', 'ACLMS')
% 
% subplot(1,3,3); hold on
% plot(real(v_h), imag(v_h), '.', 'markerSize', 3)
% plot(real(y_est_clms_h), imag(y_est_clms_h), '.', 'markerSize', 3)
% plot(real(y_est_aclms_h), imag(y_est_aclms_h), '.', 'markerSize', 3)
% plot([-10, 10], [0, 0], '--', 'color', 'red')
% plot([0, 0], [-10, 10], '--', 'color', 'red')
% hold off; xlim([-4.5, 4.5]), ylim([-4.5, 4.5])
% title('High Wind Speed', 'fontWeight', 'normal')
% legend('Data', 'CLMS', 'ACLMS')
% 
% 
