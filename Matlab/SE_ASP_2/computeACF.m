function [ b, ub ] = computeACF( sig )
%Computes biased and unbiased ACF

b = xcorr(sig, 'biased');
ub = xcorr(sig, 'unbiased');

end

