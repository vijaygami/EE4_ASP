function [ x ] = x_axis( N, F )
% Generates vector of x-axis values so they can be used when plotting FFT or ACF centred around 0.
% N = Number of samples, F = Max freq. For example, use pi when making a normalised frequency axis
    if mod(N,2) == 0   x=[-(N/2):1:(N/2)-1]*2*F/N;        % N is Even
    else               x=[-((N-1)/2):1:((N-1)/2)]*2*F/N;  % N is Odd
    end
end

