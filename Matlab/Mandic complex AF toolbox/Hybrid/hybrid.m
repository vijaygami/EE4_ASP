% FUNCTION hybrid() implement the hybrid filter combining the CLMS
% algorithm and ACLMS algorithm.
% 
% Based on the paper "Collaborative adaptive filtering in the complex domain"
% IEEE Workshop on Machine Learning for Signal Processing, 421-425, 2008.
%
% INPUT:
% x: input signal 
% N: filter length of CLMS (that of ACLMS is N1 * 2)
% mu1: step-size of CLMS
% mu2: step-size of ACLMS
% lambda: mixing parameter
% mu3: step-size of adaptation of lambda
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
function y = hybrid(x,N,mu1,mu2,lambda,mu3)

M = 1;% prediction horizon
L = length(x)-M; % 
filterinput = zeros(N,L); 
WCLMS = zeros(N,1); % weight vector of subfilter 1
eCLMS = zeros(1,L); % error of subfilter 1
suboutput1 = zeros(1,L); % output of subfilter 1

WACLMS = zeros(2*N,1); % weight vector of subfilter 2
eACLMS = zeros(1,L); % error of subfilter 2
suboutput2 = zeros(1,L); % output of subfilter 2

ehybrid = zeros(1,L);
filteroutput = zeros(1,L); 

for i = 1:L
    for m = 1:N
        if (i-m+1)>0
            filterinput(m,i) = x(1,i-m+1);
        else
            filterinput(m,i) = 0;
        end
    end
    suboutput1(i) =  transpose(filterinput(:,i)) * WCLMS;
    suboutput2(i) =  transpose([filterinput(:,i);conj(filterinput(:,i))]) * WACLMS;
    eCLMS(i) = x(1,i+M) - suboutput1(i);
    eACLMS(i) = x(1,i+M) - suboutput2(i);
    filteroutput(i) = lambda * suboutput1(i) + (1-lambda) * suboutput2(i);
    ehybrid(i) = x(1,i+M) - filteroutput(i);% error(k)
    WCLMS = WCLMS + mu1 * eCLMS(i) * conj(filterinput(:,i));% weight update of subfilter 1
    WACLMS = WACLMS + mu2 * eACLMS(i) * conj([filterinput(:,i);conj(filterinput(:,i))]);% weight update of subfilter 2
    lambda = lambda + mu3 * real(ehybrid(i) * conj(suboutput1(i)-suboutput2(i)));% update of lambda
end
y = filteroutput;

