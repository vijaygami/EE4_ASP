clc;clear;
N=2000;
fs=1000;
fo=50;
phi=0;
u=0.1;

% Generate the complex clarke voltage for different fault conditions
v(1,:)=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), 0, 0);                                    % Balanced Normal
v(2,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.464*ones(1,N), 0.464*ones(1,N), 0, 0);                  % type A
v(3,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), ones(1,N), ones(1,N), 0, 0);                              % type B
v(4,:)=clarke(N, fs, fo, phi, ones(1,N), 0.707*ones(1,N), 0.707*ones(1,N), -pi/12, pi/12);               % type C
v(5,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.896*ones(1,N), 0.896*ones(1,N), pi/12, -pi/12);         % type D
v(6,:)=clarke(N, fs, fo, phi, ones(1,N), 0.464*ones(1,N), 0.464*ones(1,N), 0, 0);                        % type E
v(7,:)=clarke(N, fs, fo, phi, 0.464*ones(1,N), 0.75*ones(1,N), 0.75*ones(1,N), pi/15, -pi/15);           % type F
v(8,:)=clarke(N, fs, fo, phi, 0.821*ones(1,N), 0.574*ones(1,N), 0.574*ones(1,N), -pi*0.0867, pi*0.0867); % type G

for i=1:8
    
    [~, h_clms, ~ ] = CLMS_AR(v(i,:), u, 1);
    [~, h_aclms, g_aclms, ~ ] = ACLMS_AR(v(i,:), u, 1);
    h_clms=conj(h_clms);
    h_aclms=conj(h_aclms);
    g_aclms=conj(g_aclms);
    f_clms(i,:)=(fs/(2*pi))*atan(imag(h_clms)./real(h_clms));
    f_aclms(i,:)=(fs/(2*pi))*atan((sqrt((imag(h_aclms).^2)-abs(g_aclms).^2))./real(h_aclms));

    % Slightly different method
    a1=(-1j.*imag(h_aclms)-1j.*sqrt(((imag(h_aclms)).^2) - (abs(g_aclms).^2)))./g_aclms;
    f_aclms1(i,:) = (fs/(2*pi))*asin(imag(h_aclms+a1.*g_aclms));
    
end


figure(1); clf; 
titles={'Frequency Estimate: No Fault'; 'Frequency Estimate: Type A Sag'; 'Frequency Estimate: Type B Sag'; 'Frequency Estimate: Type C Sag'; 'Frequency Estimate: Type D Sag'; 'Frequency Estimate: Type E Sag'; 'Frequency Estimate: Type F Sag'; 'Frequency Estimate: Type G Sag' };

for i=1:8
    subplot(2,4,i); hold on
    plot(abs(f_clms(i,:)))
    plot(abs(f_aclms(i,:)))
    box off; grid on; hold off;
    ylim([0, 100]); xlim([0, 300])
    title(titles{i}, 'fontWeight', 'normal')
    xlabel('Time [samples]')
    ylabel('Frequency [Hz]')
    legend('CLMS', 'ACLMS')
end

%% plots of how steady state frequency estimate varies with imbalance.
N=2000;
fs=1000;
fo=50;
phi=0;
u=0.1;
v=[]; 
f_aclms=[];
f_clms=[];
k=100;   % num of points to  plot
f_clms_ss=zeros(k,2);
f_aclms_ss=zeros(k,2);
delta_b=linspace(-pi*2/6, pi*2/6, k);
V_a=linspace(0, 2, k);

for i=1:k
    
    % phase imbalance
    v=clarke(N, fs, fo, phi, ones(1,N), ones(1,N), ones(1,N), delta_b(i), 0); 
    [~, h_clms, ~ ] = CLMS_AR(v, u, 1);
    [~, h_aclms, g_aclms, ~ ] = ACLMS_AR(v, u, 1);
    h_clms=conj(h_clms);
    h_aclms=conj(h_aclms);
    g_aclms=conj(g_aclms);
    f_clms=(fs/(2*pi))*atan(imag(h_clms)./real(h_clms));
    f_aclms=(fs/(2*pi))*atan((sqrt((imag(h_aclms).^2)-abs(g_aclms).^2))./real(h_aclms));
    f_clms_ss(i,1)=mean(f_clms(1750:end));      % steady state values
    f_aclms_ss(i,1)=mean(f_aclms(1750:end));

    % magnitude imbalance
    v=clarke(N, fs, fo, phi, ones(1,N).*V_a(i), ones(1,N), ones(1,N), 0, 0); 
    [~, h_clms, ~ ] = CLMS_AR(v, u, 1);
    [~, h_aclms, g_aclms, ~ ] = ACLMS_AR(v, u, 1);
    h_clms=conj(h_clms);
    h_aclms=conj(h_aclms);
    g_aclms=conj(g_aclms);
    f_clms=(fs/(2*pi))*atan(imag(h_clms)./real(h_clms));
    f_aclms=(fs/(2*pi))*atan((sqrt((imag(h_aclms).^2)-abs(g_aclms).^2))./real(h_aclms));
    f_clms_ss(i,2)=mean(f_clms(1750:end));      % steady state values
    f_aclms_ss(i,2)=mean(f_aclms(1750:end));
    
    
end
%%
figure(2); clf; 
subplot(1,2,1)
hold on;
plot(delta_b, f_clms_ss(:,1))
plot(delta_b, f_aclms_ss(:,1))
hold off
legend('CLMS', 'ACLMS', 'location', 'best')
xlabel('$\Delta_b$ [Radians]')
ylabel('Frequency Estimate [Hz]')
title('Steady State Frequency Estimate with Phase Imbalance', 'fontWeight', 'normal')
box off; grid on;
xlim([-1.1, 1.1])
ylim([35, 55])

subplot(1,2,2)
hold on;
plot(V_a, f_clms_ss(:,2))
plot(V_a, f_aclms_ss(:,2))
hold off
legend('CLMS', 'ACLMS', 'location', 'best')
xlabel('$V_a$ [Per Unit]')
ylabel('Frequency Estimate [Hz]')
title('Steady State Frequency Estimate with Magnitude Imbalance', 'fontWeight', 'normal')
box off; grid on; 

