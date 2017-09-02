% This program implements the Adaptive Amplitude CRTRL (AACRTRL).
%
% Based on the paper "Nonlinear adaptive prediction of complex-valued signals by complex-valued PRNN",
% IEEE Transactions on Signal Processing, vol 53, no 5, 2005.
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
%
%
% ...........................................

close all;
clear all;

M=4; % number of input tap (feedforward order) 
N=2; % number of neurons (feedback order)
eta=0.1; % learning rate
beta=1; % activation slope
phi=0.1; % step size
len=2000; % input signal sample size

% initialised all the variables to zeros
w=zeros(M+N+1,N); % weights of neuron
lamda=zeros(N,1); % adaptive amplitude
PI_old=zeros(M+N+1,N,N); % previous sample PI
PI_new=zeros(M+N+1,N,N); % current sample PI
Y_old=zeros(N,1); % Previous output matrix
Y_out=zeros(N,1); % output matrix
x=zeros(M,1); % input samples
% Eave=zeros(1,len);
Elog=zeros(1,len); % cost function

for monte=1:30 % monte carlo simulation
    
monte

%weight initialization
W_init=0.01*rand(M+N+1,N); %initialise the weight randomly
w=W_init;

%initial lambda (adaptive amplitude)
lamda(:,1)=1;

% generate input data
noise=wgn(1,len,1); % white gaussian noise
wgnoise=nonlinearfilt(noise); % passed the WGN into the nonlinear filter
d=0.1*(wgnoise/max(wgnoise))+0.1; % scale the desired signal

% prediction signal
xin(1)=0;
xin(2:1:len)=d(1:1:len-1); %input signal

for k=1:len
    
    x=[xin(k);x(1:M-1)]; % input sample to the network based on input tap size
    Y_old=[Y_out(1);Y_old(1:N-1)];
    Uin=[Y_old;1;x]; % input to the network 
    Vout=w'*Uin;
    
    % First Activation Function (sigmoid)
    sig_function = lamda(:,k) ./ ( 1 + exp( -( beta .* Vout )));

    % First derivative Activation Function 
    sig_function_der = beta .* sig_function .* ( 1 - sig_function );
    
    Y_out=sig_function; % value at the output of the neuron
    %Y_old=Y_out; % previous sample
    
    e(:,k)=d(:,k)-Y_out(1,:); % error 
    E(:,k)=(1/2)*(e(:,k)).^2; % cost function
    
    % Adaptive amplitude
    lamda(:,k+1)=lamda(:,k) + phi.*e(:,k).*(sig_function ./ lamda(:,k));
    
    E_dB(k)=10*log10(E(:,k)); % evaluation of cost function in dB
    
    % Compute values of pi
    
    for r=1:1+M+N
        for s=1:N
            for t=1:N
                temp=0;
                for i=1:N
                    temp=temp+w(r,i).*PI_old(i,s,t);
                end
                PI_new(r,s,t)=(sig_function_der(t,:)).*(Uin(r)+temp);
            end
        end
    end
    
 PI_old=PI_new;
 
 % Calculate weight changes
 dW=zeros(N,M+N+1);
 dW=eta.*e(:,k).*PI_new;
 w=w+dW(:,:,1);
 
end % for k

% Eave=Eave+E;
Elog=Elog+E_dB;

end % monte

Elog=Elog/monte;
% Eave_dB=10*log10(Eave/monte);

%plot graph
figure(1);
plot (1:len,d(1:len),'r',1:len,xin(1:len),'b')
legend('Target', 'Input');
figure(2);
plot(e)
title('e(k)=d(k)-y(k)')
figure(3);
plot(E)
title('E')
% axis([0 500 -0.01 0.5]);
figure(4);
plot(Elog)
title('Elog')
% figure(5)
% plot(Eave_dB)
% title('Eave')