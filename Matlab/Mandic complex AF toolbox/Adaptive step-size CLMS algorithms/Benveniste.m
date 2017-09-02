% FUNCTION Benveniste() implements the Benveniste's variable step size CLMS algorithm 
%
% Based on the book of A. Benveniste, M. Metivier and P. Priouret, Adaptive Algorithms and Stochastic Approximation, New York: Spinger-Verlag, 1990
%
% INPUT:
% x: input signal which should be scaled according to the dynamic range of nonlinearity 
% N: filter length
% mu: step-size
% rho: step-size of adaptation of mu
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
function y = Benveniste(x,N,mu,rho)

M = 1;% prediction horizon
L = length(x)-M; % length of simulation
filterinput = zeros(N,L);%input of FIR
filteroutput = zeros(1,L);%output of FIR
learning = zeros(1,L);% the adaptive learning rate;
Phi = zeros(N,1);
WVSLMS = zeros(N,1);% weight
eVSLMS = zeros(1,L);% error
EVSLMS = zeros(1,L);% mean square error
filteroutput = zeros(1,L);% output


for i = 1:L
    for m = 1:N
        if (i-m+1)>0
            filterinput(m,i) = x(1,i-m+1);
        else
            filterinput(m,i) = 0;
        end
    end 

    filteroutput(i) = transpose(filterinput(:,i)) * WVSLMS;%
    eVSLMS(i) = x(i+M) - filteroutput(i);% 
    EVSLMS(i) = 10 * log10(1/2 * eVSLMS(i)' * eVSLMS(i));
    if i == 1
        learning(1) = mu;
        Phi = zeros(N,1);
    else
        Phi =  (eye(N,N) - learning(i-1) * conj(filterinput(:,i-1)) * transpose(filterinput(:,i-1))) * Phi + ...
                eVSLMS(i-1) * conj(filterinput(:,i-1));
        learning(i) = learning(i-1) + rho * real(eVSLMS(i) * filterinput(:,i)'* conj(Phi));
    end
    WVSLMS = WVSLMS + learning(i) * eVSLMS(i) * conj(filterinput(:,i));
end
y = filteroutput;



