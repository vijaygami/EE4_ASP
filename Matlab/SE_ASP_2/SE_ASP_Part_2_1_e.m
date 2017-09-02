clc;clear;clf;
rng(1);

L=256; 

% MUSIC method, N=30
N=30; n=0:N-1;     
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;
[X,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,2,[],1,'corr');

subplot(1,3,1);
plot(F,S); set(gca,'xlim',[0.25 0.40]);
grid on; box off;  xlabel('Frequency [Hz]'); ylabel('Pseudospectrum');
title('MUSIC Method, N = 30', 'fontWeight', 'normal')
ylim([0, 750])


%Music method N=25
rng(1)
N=25;       
n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;
[X,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,2,[],1,'corr');

subplot(1,3,2);
plot(F,S); set(gca,'xlim',[0.25 0.40]);
grid on; box off; xlabel('Frequency [Hz]'); ylabel('Pseudospectrum');
title('MUSIC Method, N = 25', 'fontWeight', 'normal')
ylim([0, 500])


% compare with periodgram method, N=30
rng(1);
N=30; n=0:N-1;
noise=0.2/sqrt(2)*(randn(1,N)+1j*randn(1,N));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;
subplot(1,3,3);
plot(x_axis(L, 0.5)+0.5, (fftshift(pgm(x, L))))
xlabel('Frequency [Hz]'); ylabel('PSD [Power/Hz]');
grid on; box off; xlim([0.25 0.40]); ylim([0, 40])
title('Periodogram, N = 30', 'fontWeight', 'normal')


