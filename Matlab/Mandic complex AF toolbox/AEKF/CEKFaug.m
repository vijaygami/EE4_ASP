% FUNCTION CEKFaug performs the complex-valued EKF algorithm.
%
% Based on the paper "An augmented Extended Kalman Filter algorithm for complex-valued Recurrent Neural Networks", Neural Computation, vol 19, no 4, 2007.
%
% INPUT:
% input: input signal which should be scaled according to the dynamic range of nonlinearity 
% target: desired response
% mode: choice of nonlinearity
% beta: the value of slope of nonlinearity
% p: length of input vector
% KalmanR: covariance matrix
% KalmanQ: covariance matrix
% KalmanP: covariance matrix
% monte: number of iterations

% OUTPUT:
% y: output signal
% e: output error signal
% w: weight vector of adaptive filter
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
function [y, e, W, weightRecord_aug,K,A,OutputVar,P,W_aug] = CEKFaug(input, target, beta, mode, p, N, KalmanR, KalmanQ, KalmanP, monte)

pp = length(input);            % Number of input signal
target = target + KalmanR*(randn(size(target))+i*randn(size(target))); % network target vector

s = zeros(p,1) + i*zeros(p,1);                         % The input vector 
DW = zeros(p+N+1, N) + i*zeros(p+N+1, N);                    % Change applied to the weight matrix for neurons 1&2



% Weights vectors
W = 0.01*randn(p+N+1,N) + i*0.01*randn(p+N+1,N);                     % The weight matrix
W_conj= conj(W);
W_aug = [W;W_conj];
PI_prev_aug = zeros(2*(p+N+1),N*N) + i*zeros(2*(p+N+1),N*N);  PI_curr_aug = zeros(2*(p+N+1),N*N) + i*zeros(2*(p+N+1),N*N);  % Jacobian matrix

K = zeros(N*N, 2*(p+N+1))+ i*zeros(N*N, 2*(p+N+1));
P = KalmanP*(eye(2*(p+N+1),2*(p+N+1))+i*eye(2*(p+N+1),2*(p+N+1)));
R = KalmanR*(eye(N*N)+i*eye(N*N));
Q = KalmanQ*(eye(2*(p+N+1),2*(p+N+1))+i*eye(2*(p+N+1),2*(p+N+1)));

weightRecord_aug = zeros(2*(p+N+1),pp) + i*zeros(2*(p+N+1),pp);


E = zeros(pp,1) + i*zeros(pp,1);
MSE = zeros(pp,1) + i*zeros(pp,1);
% Output vectors from the neurons
Y_curr_aug = zeros(2*pp, N) + i*zeros(2*pp, N);

U = zeros(p+N+1, 1) + i*zeros(p+N+1, 1);              % INput to the network(external,feedback and bias) 
U_conj = conj(U);
V = zeros(1,N) + i*zeros(1,N);                        % weighted input vector to a neuron;
Errors = 0.001*(randn(pp, 1) + i*randn(pp, 1));
Phi = zeros(1, N) + i*zeros(1,N);

for Monte=1:monte
for k=2:pp
    
    s = [input(k-1);s(1:p-1)];
    U(:) = [s;1+i;Y_curr_aug(k-1,:).'];      % The input vector to the network
    U_conj(:) = conj(U(:));
    U_aug = [U(:);U_conj(:)];

   V_aug = U_aug.'*W_aug;
   % keyboard;
    switch (mode)
        case 1
            %Y_curr(k,:)     = logsig(beta*V);
             Y_curr_aug(k,:)    = 1./(1+exp(-beta*V_aug));
             
        case 2
            %Y_curr(k,:)     = tansig( beta * V);
             Y_curr_aug(k,:)     = tanh( beta * V_aug);
             %Y_curr     = tanh( beta * V)
    end

    Errors(k)     = target(k)-Y_curr_aug(k,1) ;   %taking the delay of observations to minus the y output
    
	switch(mode)
        case 1
        %Phi = dlogsig(beta*V,Y_curr(k,:));
        Phi = beta * (1./(1+exp(-beta*V_aug))).*( 1 - 1./(1+exp(-beta*V_aug)) );
        case 2
        %Phi = dtansig(beta * V,Y_curr(k,:)); 
        Phi = 1 - (tanh(beta * V_aug)).^2;
         %Phi_check = 1 - (tanh(beta * V)).^2
	end

	for j=1:N      %Put N PI matrix into one Big PI, that's N x [(1+p+N) x N] size
        for n=1:N
            for l=(1:2*(p+1+N))
                temp = 0;
                for m=1:N
                    temp = temp + conj(W_aug(1+p+m, j))*PI_prev_aug(l,(m-1)*N+n);
                end
                wtemp(n,l) = temp;
                if n==j
                    wtemp(n, l) = wtemp(n, l) + conj(U_aug(l));
                end
                PI_curr_aug(l,(j-1)*N+n) = conj(Phi(1,j)) * wtemp(n, l);
            end
        end
    end


    PI_prev_aug = PI_curr_aug;
    
    A = (R+conj(PI_curr_aug.')*(P+Q)*(PI_curr_aug)); % danger of matrix singularity due to inversion

    K = (P+Q)*(PI_curr_aug)*(A^(-1));
    W_aug = W_aug + K(:,1:N:N*N)*Errors(k);
    P = P - K*conj(PI_curr_aug.')*(P+Q)+Q;
    
    
    OutputVar(:,k) = diag(R+conj(PI_curr_aug.')*(P)*(PI_curr_aug));
    P_record(:,:,k) = P;
    weightRecord_aug(:, k) = W_aug(:,1);
end

E = E +(1/2)*abs(Errors).^2;
MSE = MSE+10*log10((1/2)*abs(Errors).^2);
var_error(Monte) = var(Errors(1:end));
var_signal(Monte) = var(input(1:end));
Rp(Monte) = 10*log10(var_signal(Monte)/var_error(Monte));

end
%keyboard;
y = Y_curr_aug(1:pp,1);

e = Errors;

  error = norm(y(1:end).' - input(1:end)); % check on this

  fprintf('EKF error    = %6.2f',error)
  fprintf('\n')
  
  fprintf('Rp     = %6.2f',mean(Rp))
  fprintf('\n')
  
  % PLOT RESULTS:
  % ============
											

 figure(1)

  plot(1:length(target),abs(target),'k+:',1:length(y),abs(y),'r',1:length(input),abs(input),'b');
  ylabel('Prediction','fontsize',15)
  %axis([0 pp -1.5 1.5]); 

  figure(2)
  plot(1:length(OutputVar),mean(abs(OutputVar)),'m');
  ylabel('Output variance','fontsize',15);
  %axis([0 pp 1.4 2]);
  
  figure(3)
  plot(E)
  
  figure(4)
  plot(MSE)

