function [ y, fn ] = fm_gen( fs, noise_var)
% Generates the FM signal given by equation 50 page 18 CW handout

    rng(1);

    n=[501:1000];
    n2=[1001:1500];
    fn = [100*ones(1,500), 100*ones(1,500)+(n-500)/2, 100*ones(1,500)+(((n2-1000)/25).^2)];
    phi=cumsum(fn);

    eta = wgn(1, 1500, 10*log10(noise_var),'complex');   % Additive complex noise
    y = exp(1i*(2*pi*phi./fs)) + eta;

end

