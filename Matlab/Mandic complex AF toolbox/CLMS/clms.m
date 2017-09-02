% FUNCTION clms() implements the CLMS algorithm 
% 
% Based on Table 5.1, Adaptive Filter Theory, Simon Haykin.
%
% INPUT:
% x: input signal 
% N: filter length
% mu: step-size
%
% OUTPUT:
% y: filter output

% ...........................................
function y = clms(x,N,mu)

M = 1;% prediction horizon
L = length(x)-M;
filterinput = zeros(N,L); 
WLMS = zeros(N,1); % weight 
eLMS = zeros(1,L); % error
filteroutput = zeros(1,L); 

for i = 1:L
    for m = 1:N
        if (i-m+1)>0
            filterinput(m,i) = x(1,i-m+1);
        else
            filterinput(m,i) = 0;
        end
    end
    filteroutput(i) = transpose(filterinput(:,i)) * WLMS;% filter output
    eLMS(i) = x(1,i+M) - filteroutput(i);% error(k)
    WLMS = WLMS + mu * eLMS(i) * conj(filterinput(:,i));% weight update
end
y = filteroutput;
