% Delay Vector Variance method for real and complex signals
%
%
% USAGE: C = dvv (X, m, Nsub, nd, Ntv)
%
% Input:
%    X       original real-valued or complex time series
%    m       delay embedding dimension
%    Ntv     number of points on horizontal axes
%    Nsub    number of reference DVs to consider
%    nd      Span over which to perform DVV
% ...........................................
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
function data = dvv (X, m, Nsub, nd, Ntv)

% Default parameters
if (nargin<1)
    error('Not enough Input arguments');
end
if (nargin<2)
    m = 4;
end
if (nargin<3)
    Nsub = 200;
end
if (nargin<4)
    nd = 2.0;
end
if (nargin<5)
    Ntv = 25*nd;
end

% Initial Conditions
N = length(X);              % Length of input vector
tau = 1;                    % Time delay parameter
d = zeros(N-m*tau, Nsub);
y = zeros(Ntv,1);
count = 0;
acc = 0;


% Makes input vector X a column vector
if (size(X,2)>size(X,1))
    X = X';
end

% Generate Nsub subset from existing DV's, randomly
temp = randperm (N - m*tau);
ref = temp(1:Nsub) + m*tau;

% Computes the pair wise distances b/w reference DV's and all DV's
for i = 1:Nsub
    for j = m*tau+1:N
        d(j-m*tau,i) = norm (X(ref(i)-m*tau:tau:ref(i)-tau) - X(j-m*tau:tau:j-tau));
        if (ref(i) ~= j)
            acc = acc + d(j-m*tau,i);
            count = count + 1;
        end
    end
end

% Mean and std variation calculation of input data
avg = acc/count;
count = 0;
acc = 0;
for i = 1:Nsub
    for j = m*tau + 1:N
        if (ref(i) ~= j)
            acc = acc + (d(j-m*tau,i)-avg).^2;
        end
    end
end
variance = sqrt(acc/(count-1));

% Calculates the range vector consisting of Ntv equally spaced regions
n = (1:Ntv)-1;
rd = avg-nd*variance + (2*nd*variance*n)/(Ntv-1);

% Creates sets of DV's, for each ref element of subset and value rd, which have norms closer than distance rd to ref
for n = 1:length(rd)
    
    if (rd(n)>0)
        
        tot = 0;
        count = 0;
        
        for k=1:Nsub
            
            IND = find(d(:,k) <= rd(n)) + m*tau;
            IND = IND(IND~=k);
            
            
            % Only those variance values are considered for which the corresponding
            % sets have atleast 30 DVs
            if (length(IND) >= 30)
                tot = tot + var(X(IND));
                count = count+1;
            end
        end
        if (~count)
            y(n) = NaN;
        else
            y(n) = tot/(count*var(X));
        end
        
    else
        y(n) = NaN;
    end
    
end

% Horizontal axis
T = (rd'-avg)/variance;

% DVV Output
data = [T y];
