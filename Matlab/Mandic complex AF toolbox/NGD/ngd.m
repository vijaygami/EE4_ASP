% FUNCTION ngd() implements the operation of NGD algorithm 
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
function y = ngd(x,N,mu)

M = 1;% prediction horizon
L = length(x)-M; % run length
filterinput = zeros(N,L); 
WNGD = zeros(N,1); % weight vector 
eNGD = zeros(1,L); % error
filteroutput = zeros(1,L); 
output = zeros(1,L);

for i = 1:L
    for m = 1:N
        if (i-m+1)>0
            filterinput(m,i) = x(1,i-m+1);
        else
            filterinput(m,i) = 0;
        end
    end % inputing FIR
    filteroutput(i) = transpose(filterinput(:,i)) * WNGD;% output of FIR filter
    output(i) = f(filteroutput(i));
    eNGD(i) = x(i+M) - output(i);% error(k) of nonlinear FIR filter
    WNGD = WNGD + mu * eNGD(i) * conj(fderive(filteroutput(i))) * conj(filterinput(:,i));% weight update
end
y = output; 


