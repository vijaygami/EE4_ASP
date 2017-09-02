% FUNCTION cekf performs the complex-valued EKF algorithm.
%
% Based on the paper of SuLee, "An augmented Extended Kalman Filter algorithm for complex-valued Recurrent Neural Networks", Neural Computation, vol 19, no 4, 2007.
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
%
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
function [y, e, W, weightRecord,K,A,OutputVar,P] = cekf(input, target, beta, mode, p, N, KalmanR, KalmanQ, KalmanP, monte)

pp = length(input);            % Number of input signal

%ita=0.1
target = target + KalmanR*(randn(size(target))+i*randn(size(target))); % network target vector

s = zeros(p,1) + i*zeros(p,1);                         % The input vector 
%DW = zeros(p+N+1, N) + i*zeros(p+N+1, N);                    % Change applied to the weight matrix for neurons 1&2

% Weights vectors
W = eye(p+N+1,N) + i*eye(p+N+1,N);                     % The weight matrix
PI_prev = zeros(N*N, p+N+1) + i*zeros(N*N, p+N+1);  PI_curr = zeros(N*N, p+N+1)+ i*zeros(N*N, p+N+1);  % Jacobian matrix
K = zeros(N*N, p+N+1)+ i*zeros(N*N, p+N+1);
P = KalmanP*(eye(p+N+1,p+N+1)+i*eye(p+N+1,p+N+1));
R = KalmanR*(eye(N*N)+i*eye(N*N));
Q = KalmanQ*(eye(p+N+1,p+N+1)+i*eye(p+N+1,p+N+1));

weightRecord = zeros(p+N+1,pp) + i*zeros(p+N+1,pp);
P_record = zeros(p+N+1,p+N+1,pp) + i*zeros(p+N+1,p+N+1,pp);
OutpurVar = zeros(1,pp) + i*zeros(1,pp);
E = zeros(pp,1) + i*zeros(pp,1);
MSE = zeros(pp,1) + i*zeros(pp,1);
% Output vectors from the neurons
Y_curr = zeros(pp, N) + i*zeros(pp, N);

U = zeros(p+N+1, 1) + i*zeros(p+N+1, 1);              % INput to the network(external,feedback and bias)     
V = zeros(1,N) + i*zeros(1,N);                        % weighted input vector to a neuron
Errors = 0.001*(randn(pp, 1) + i*randn(pp, 1));
Phi = zeros(1, N) + i*zeros(1,N);

for Monte=1:monte
for k=2:pp
    s = [input(k-1);s(1:p-1)];
    U(:) = [s;1;Y_curr(k-1,:).'];      % The input vector to the network
    V  = U.'*W;
        
    switch (mode)
        case 1
            %Y_curr(k,:)     = logsig(beta*V);
             Y_curr(k,:)     = 1./(1+exp(-beta*V));
        case 2
            %Y_curr(k,:)     = tansig( beta * V);
             Y_curr(k,:)     = tanh( beta * V);
    end

    Errors(k)     = target(k)-Y_curr(k, 1) ;   %taking the delay of observations to minus the y output
    Errors_ori(k)   = (input(k)-Y_curr(k, 1)) ;   %taking the delay of observations to minus the y output

	switch(mode)
        case 1
        %Phi = dlogsig(beta*V,Y_curr(k,:));
        Phi = beta * (1./(1+exp(-beta*V))).*( 1 - 1./(1+exp(-beta*V)) );
        case 2
        %Phi = dtansig(beta * V,Y_curr(k,:)); 
        Phi = 1 - (tanh(beta * V)).^2;
	end

	for j=1:N      %Put N PI matrix into one Big PI, that's N x [(1+p+N) x N] size
        for n=1:N
            for l=1:p+1+N
                temp = 0;
                for m=1:N
                    temp = temp + conj(W(1+p+m, j))*PI_prev((m-1)*N+n, l);
                end
                wtemp(n,l) = temp;
                if n==j
                    wtemp(n, l) = wtemp(n, l) + conj(U(l));
                end
                PI_curr((j-1)*N+n, l) = conj(Phi(j)) * wtemp(n, l);
            end
        end
    end
    PI_prev = (PI_curr);

    A = (R+conj(PI_curr)*(P+Q)*(PI_curr.')); % danger of matrix singularity due to inversion

    K = (P+Q)*(PI_curr.')*(A^(-1));
    W = W + K(:,1:N)*Errors(k);
    P = P - K*conj(PI_curr)*(P+Q)+Q;
%keyboard;
    OutputVar(:,k) = diag(R+conj(PI_curr)*(P)*(PI_curr.'));
    weightRecord(:, k) = W(:,1);
    P_record(:,:,k) = P;
end

E = E +(1/2)*abs(Errors).^2;
MSE = MSE+10*log10((1/2)*abs(Errors).^2);
var_error(Monte) = var(Errors(1:end));
var_signal(Monte) = var(input(1:end));
Rp(Monte) = 10*log10(var_signal(Monte)/var_error(Monte));

end
%keyboard;
y = Y_curr(:,1);
e = Errors;

  errorekf = norm(y(1:end).' - input(1:end)); % check on this

  fprintf('EKF error    = %6.2f',errorekf)
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
  save Ewind_ekf E
  
  figure(4)
  plot(MSE)
