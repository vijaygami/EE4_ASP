clc;clear; clf;

L=256;

%_ACF 'x' as given by equation 12 in section 1.2 of courcework_____________________________________
M=10;
x1=(1/M)*[[M:-1:1],zeros(1,(L-(2*M-1))),[1:1:M-1]];  % Zero padded vector of length L

M=128;
x2=(1/M)*[[M:-1:1],zeros(1,(L-(2*M-1))),[1:1:M-1]];  % Zero padded vector of length L

%_Alternate forms of 'x' to jutify the original choice____________________________________________
M=10;
x3=(1/M)*[[M:-1:1], [1:1:M-1],zeros(1,(L-(2*M-1)))];  % No circular shifting when padding

x4=(1/M)*[[1:1:M-1],[M:-1:1], zeros(1,(L-(2*M-1)))];  % Linear dhifting

x5=[[M:-1:1],0, [1:1:M-1]]; % No padding and length 20



subplot(2, 5, 1)
plot(x1)
title({'Circular Shifted ACF', 
    '~~~~~M=10, L=256'},'FontWeight', 'normal')

subplot(2, 5, 6)
plot(x_axis(L, 0.5) , fftshift((fft(x1))))
title('FFT(ACF)', 'FontWeight', 'normal')

subplot(2, 5, 2)
plot(x2)
title({'Circular Shifted ACF', 
   '~~~~M=128, L=256'}, 'FontWeight', 'normal')

subplot(2, 5, 7)
plot(x_axis(L, 0.5) , fftshift((fft(x2))))
title('Real(FFT(ACF))', 'FontWeight', 'normal')

subplot(2, 5, 3)
plot(x3)
title({'Non-Symetric ACF', 
    '~~~~M=10, L=256'}, 'FontWeight', 'normal')

subplot(2, 5, 8)
plot(x_axis(L, 0.5) , fftshift((fft(x3))))
title('Real(FFT(ACF))', 'FontWeight', 'normal')

subplot(2, 5, 4)
plot(x4)
title({'Linear Shifted ACF', 
    '~~~~M=10, L=256'}, 'FontWeight', 'normal')

subplot(2, 5, 9)
plot(x_axis(L, 0.5) , fftshift((fft(x4))))
title('Real(FFT(ACF))', 'FontWeight', 'normal')

subplot(2, 5, 5)
plot(x5)
title({'Circular shifted ACF', 
    '~~~~~M=10, L=20'}, 'FontWeight', 'normal')

subplot(2, 5, 10)
plot(x_axis(20, 0.5) , fftshift((fft(x5))))
title('FFT(ACF)', 'FontWeight', 'normal')


for i=1:5
    subplot(2, 5, i)
    box off; 
    grid on; 
    xlabel('Sample Number')
    ylabel('ACF')
end

for i=6:10
    subplot(2, 5, i)
    box off; 
    grid on; 
    ylabel('FFT(ACF)')
    xlabel({'Normalised Frequency' 
    '~~~~~[cycles/sample]'})
end
