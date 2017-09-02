% periodogram function withount using fft function
% Modified from my year 3 2015 ASP courcework

function [pxx]= pgm(x, L);
  %L = length of FFT
  
  N=length(x);
  f=fftshift(fft(x, L));
  pxx=(f.*conj(f))./N;
end



