function [ x ] = sine_gen( freq, fs, N )

t=[0:(N-1)];
x = sin(2*pi*freq*t/fs);

end
