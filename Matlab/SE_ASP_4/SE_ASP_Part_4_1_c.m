clc;clear;
N=200;
fs=10000;
fo=50;
phi=0;

s=3; % number of subplots

% Unbalanced with different phases
figure(1); clf;
subplot(1,s,1); hold on
v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), pi/6, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, -pi/6);
plot(real(v), imag(v), '.', 'markerSize', 5)
hold off;
legend('Balanced','$\Delta_b$=$\pi$/6, $\Delta_c$=0', '$\Delta_b$=0, $\Delta_c$=-$\pi$/6', 'location', 'south')
title('Effects of Phase Imbalance', 'fontWeight', 'normal')


% Unbalanced with different mangnitude
subplot(1,s,2); hold on
v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N)*0.5, ones(1,N)*1, ones(1,N)*1, 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N)*1, ones(1,N)*0.5, ones(1,N)*1, 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N)*1, ones(1,N)*1, ones(1,N)*0.5, 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
hold off; 
title('Effects of Amplitude Imbalance', 'fontWeight', 'normal')
legend('Balanced','$v_a$=0.5, $v_b$=1, $v_c$=1', '$v_a$=1, $v_b$=0.5, $v_c$=1', '$v_a$=1, $v_b$=1, $v_c$=0.5', 'location', 'south')

% Unbalanced both phase and magnitude
subplot(1,s,3); hold on;
v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, 0);
plot(real(v), imag(v), '.', 'markerSize', 5)
v=clarke(N, fs, fo, phi, ones(1,N)*0.5+0.1*randn(1,N), 1*ones(1,N)+0.1*randn(1,N), ones(1,N)*1.5, 1, 2);
plot(real(v), imag(v), '.', 'markerSize', 5)
hold off; 
title('Time Varying & Unbalanced Amplitude', 'fontWeight', 'normal')
legend('Balanced', 'Time Varying Amplitudes', 'location', 'south')

for i=1:s
    subplot(1,s,i); hold on
    plot([-10, 10], [0, 0], '--', 'color', 'red')
    plot([0, 0], [-10, 10], '--', 'color', 'red')
    hold off; grid on; box off
    xlabel('Real Part')
    ylabel('Imaginary Part')
    xlim([-1.75, 1.75]), ylim([-2.45, 1.5])
end



%% Different Sag types for different fault conditions

vsag(1,:)=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, 0);                                    % Balanced Normal
vsag(2,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.464*ones(1,N), 0.464*ones(1,N), 0, 0);                  % type A
vsag(3,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), ones(1,N), ones(1,N), 0, 0);                              % type B
vsag(4,:)=clarke(N, fs, fo, phi, ones(1,N), 0.707*ones(1,N), 0.707*ones(1,N), -pi/12, pi/12);               % type C
vsag(5,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.896*ones(1,N), 0.896*ones(1,N), pi/12, -pi/12);         % type D
vsag(6,:)=clarke(N, fs, fo, phi, ones(1,N), 0.464*ones(1,N), 0.464*ones(1,N), 0, 0);                        % type E
vsag(7,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.75*ones(1,N), 0.75*ones(1,N), pi/15, -pi/15);           % type F
vsag(8,:)=clarke(N, fs, fo, phi, 0.821*ones(1,N), 0.574*ones(1,N), 0.574*ones(1,N), -pi*0.0867, pi*0.0867); % type G


figure(2); clf


subplot(1,3,1); hold on;
plot(real(vsag(1,:)), imag(vsag(1,:)), '.', 'markerSize', 5)
plot(real(vsag(2,:)), imag(vsag(2,:)), '.', 'markerSize', 5)
plot(real(vsag(3,:)), imag(vsag(3,:)), '.', 'markerSize', 5)
hold off;
legend('No Fault', 'Sag Type A', 'Sag Type B', 'location', 'south')


subplot(1,3,2); hold on;
plot(real(vsag(1,:)), imag(vsag(1,:)), '.', 'markerSize', 5)
plot(real(vsag(4,:)), imag(vsag(4,:)), '.', 'markerSize', 5)
plot(real(vsag(5,:)), imag(vsag(5,:)), '.', 'markerSize', 5)
hold off;
legend('No Fault', 'Sag Type C', 'Sag Type D', 'location', 'south')



subplot(1,3,3); hold on;
plot(real(vsag(1,:)), imag(vsag(1,:)), '.', 'markerSize', 5)
plot(real(vsag(6,:)), imag(vsag(6,:)), '.', 'markerSize', 5)
plot(real(vsag(7,:)), imag(vsag(7,:)), '.', 'markerSize', 5)
plot(real(vsag(8,:)), imag(vsag(8,:)), '.', 'markerSize', 5)
hold off;
legend('No Fault', 'Sag Type E', 'Sag Type F', 'Sag Type G', 'location', 'south')

for i=1:3
    subplot(1,3,i); hold on
    plot([-10, 10], [0, 0], '--', 'color', 'red')
    plot([0, 0], [-10, 10], '--', 'color', 'red')
    hold off; grid on; box off
    title('Circularity Diagram for Different Sag Types', 'fontWeight', 'normal')
    xlabel('Real Part')
    ylabel('Imaginary Part')
    xlim([-1.75, 1.75]), ylim([-2.3, 1.5])
end
