function [ v ] = clarke( N, fs, fo, phi, va, vb, vc, deltab, deltac )
% Clarke transform 

if length(va) ~= N
   error('va is not of length N');
elseif length(vb) ~= N
   error('vb is not of length N');
elseif length(vc) ~= N
   error('vc is not of length N');
end

n=[1:N];
a(n)=(sqrt(6)/6)*(va(n)+vb(n)*exp(1j*deltab)+vc(n)*exp(1j*deltac));
b(n)=(sqrt(6)/6)*(va(n)+vb(n)*exp(-1j*(deltab+2*pi/3))+vc(n)*exp(-1j*(deltac-2*pi/3)));
v(n)=a(n).*exp(1j*((2*pi*fo*n/fs)+phi))+b(n).*exp(-1j*((2*pi*fo*n/fs)+phi));

end

