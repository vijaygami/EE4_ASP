% FUNCTION scngd() implements the Split Complex NGD algorithm
%
% Based on the book "Signal Processing Techniques for Knowledge Extraction and Information Fusion"
% Springer US, 2008.
%
%
% INPUT:
% x: input signal which should be scaled according to the dynamic range of nonlinearity 
% N: filter length
% mu: step-size
%
% OUTPUT:
% y: filter output
%
%
% Complex Valued Nonlinear Adaptive Filtering toolbox for MATLAB
% Supplementary to the book:
% 
% "Complex Valued Nonlinear Adaptive Filters: Noncircularity, Widely Linear and Neural Models"
% by Danilo P. Mandic and Vanessa Su Lee Goh
% 
% (c) Copyright Danilo P. Mandic 2009
% http://www.commsp.ee.ic.ac.uk/~mandic
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
% 
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
% 
%    You can obtain a copy of the GNU General Public License from
%    http://www.gnu.org/copyleft/gpl.html or by writing to
%    Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ...........................................
function y = scngd(x,N,mu)

M = 1;% prediction horizon
L = length(x)-M; % run length
filterinput = zeros(N,L); 
RWNGD = zeros(N,1); % weight vector of standard input vector
IWNGD = zeros(N,1);
realoutput = zeros(1,L);
imagoutput = zeros(1,L);
output = zeros(1,L);
eReal = zeros(1,L);
eImag = zeros(1,L);


for i = 1:L
    for m = 1:N
        if (i-m+1)>0
            filterinput(m,i) = x(i-m+1);
        else
            filterinput(m,i) = 0;
        end
    end % inputing FIR
    realinput = real(filterinput);
    imaginput = imag(filterinput);
    realoutput(i) = transpose(realinput(:,i)) * RWNGD;
    imagoutput(i) = transpose(imaginput(:,i)) * IWNGD;
    output(i) = realoutput(i) + j * imagoutput(i);
    eReal(i) = real(x(i+M)) - f(realoutput(i));
    eImag(i) = imag(x(i+M)) - f(imagoutput(i));
    RWNGD = RWNGD + mu * eReal(i) * fderive(realoutput(i)) * realinput(:,i);
    IWNGD = IWNGD + mu * eImag(i) * fderive(imagoutput(i)) * imaginput(:,i);
end
y = output;


