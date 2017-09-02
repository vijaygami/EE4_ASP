% FUNCTION augment_complex_rtrl_A() implements the Augmented CRTRL algorithm
%
% Based on the paper "A augmented CRTRL for complex-valued recurrent neural networks", 
% Neural Networks, vol 20, issue 10, 2007.
%
% INPUT: input signal which should be scaled according to the dynamic range of nonlinearity 
%
% OUTPUT:
% Y_out1A: output signal
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
function Y_out1A =augment_complex_rtrl_A(input)
%Defining parameters
p=5; %input taps
N=3;%number of neurons
alpha=0.1;%learning rate
len=5000; %length of input sample
%Initialization of weight parameter
w1=zeros(2*(p+N+1),1)+i*zeros(2*(p+N+1),1);
w2=zeros(2*(p+N+1),1)+i*zeros(2*(p+N+1),1);
w3=zeros(2*(p+N+1),1)+i*zeros(2*(p+N+1),1);
W=zeros(2*(p+N+1),N)+i*zeros(2*(p+N+1),N);

%Initialization of weight change
dW1=zeros(1,2*(p+N+1))+i*zeros(1,2*(p+N+1));
dW2=zeros(1,2*(p+N+1))+i*zeros(1,2*(p+N+1));
dW3=zeros(1,2*(p+N+1))+i*zeros(1,2*(p+N+1));
DW=zeros(N,2*(p+N+1))+i*zeros(N,2*(p+N+1));

%Output Matrix
Y_old1=zeros(len+1,1)+i*zeros(len+1,1);
Y_out1=zeros(len,1)+i*zeros(len,1);
Y_old2=zeros(len+1,1)+i*zeros(len+1,1);
Y_out2=zeros(len,1)+i*zeros(len,1);
Y_old3=zeros(len+1,1)+i*zeros(len+1,1);
Y_out3=zeros(len,1)+i*zeros(len,1);

%Preparing the basic format fo the input signal
x=zeros(p,1);
ElogA=zeros(1,len);
E_dB=zeros(1,len);  
for monte=1:100
    monte;
     
%PI for split case 
PR_old=zeros(2*(p+N+1),N,N);%previous sample real
PR_new=zeros(2*(p+N+1),N,N);%current sample real

%weight initialization
W_init1real=0.01*rand(2*(p+N+1),1);%initialise the weight1 randomly
W_init1imag=0.01*rand(2*(p+N+1),1);
w1=W_init1real+i*W_init1imag;
W_init2real=0.01*rand(2*(p+N+1),1); %initialise the weight2 randomly
W_init2imag=0.01*rand(2*(p+N+1),1);
w2=W_init2real+i*W_init2imag;
W_init3real=0.01*rand(2*(p+N+1),1); %initialise the weight2 randomly
W_init3imag=0.01*rand(2*(p+N+1),1);
w3=W_init3real+i*W_init3imag;


%Load Data of Complex Colored Input
d=input(1:len);
xin(1)=0;
xin(2:len)=d(1:len-1);

%Activity of the Neurons
for k=1:len
    x=[xin(k);x(1:p-1)]; %input of the system due to input 
    Uina=[Y_old1(k);Y_old2(k);Y_old3(k);1+i;x]; %the main input to the system, 1 represents the bias
    Uin=[Uina;conj(Uina)];

    %First Neuron Activity
    Vout1=w1.'*Uin;
    sig_function1 = tanh(Vout1);%Output of the neuron 1 real
    sig_function_der1 = ( 1 - sig_function1^2 );%the derivative,f'
    Y_out1(k)=sig_function1;%store in the output matrix of 1st neuron
    Y_old1(k+1)=Y_out1(k);
    
   
    %Second Neuron Activity
    Vout2=w2.'*Uin;
   sig_function2 = tanh(Vout2);%Output of the neuron 1 real
    sig_function_der2 = ( 1 - sig_function2^2 );%the derivative,f'
    Y_out2(k)=sig_function2;%store in the output matrix of 1st neuron
    Y_old2(k+1)=Y_out2(k);
    
     %Third Neuron Activity
    Vout3=w3.'*Uin;
   sig_function3 = tanh(Vout3);%Output of the neuron 1 real
    sig_function_der3 = ( 1 - sig_function3^2 );%the derivative,f'
    Y_out3(k)=sig_function3;%store in the output matrix of 1st neuron
    Y_old3(k+1)=Y_out3(k);
    
    %Real and Imaginary Part of Output
    u1(k)=real(Y_out1(k));
    v1(k)=imag(Y_out1(k));
    
    %Error Calculation
    e(k) = d(k) - Y_out1(k);
    
    %error component
    e_real(k)=real(d(k))-u1(k);
    e_imag(k)=imag(d(k))-v1(k);
   
    %MSE Error Calculation
    E(k)=(1/2)*(e_real(k).^2+e_imag(k).^2);
    E_dB(k)=10*log10(E(k));%error value at k step
    
    %Matrix ready before calculating Pij
    sig_function_der=[sig_function_der1,sig_function_der2,sig_function_der3];
    
     %Calculating Pij
    for l=1:2*(N+p+1)%row
        for t=1:N%
            for j=1:N%page
                tempR=0;
                m=0;
                for r=1:2*N%rotation
                tempR=tempR+conj((W(r+m*(p+1),j))).*(PR_old(l,r-m*N,j));
                if r==N
                    m=1;
                else
                end  
                end
                if t==j
                    tempR=tempR+conj(Uin(l));
                else
                end
                PR_new(l,t,j)=conj(sig_function_der(j)).*tempR;
            end
        end
    end
    PR_old=PR_new;
    
    %weight change
    DW=alpha.*(e(k).*PR_new);
    w1=w1+DW(:,1,1);
    w2=w2+DW(:,1,2);
    w3=w3+DW(:,1,3);
    W=[w1,w2,w3];
end%k
Y_out1A=Y_out1;
ElogA=ElogA+E_dB;
var_error(monte)=var(abs(e(1000:end)));
var_signal(monte)=var(d(1000:end));
Rp(monte)=10*log10(var_signal(monte)/var_error(monte));

end%monte
ElogA=ElogA/monte;
PredictionGain=mean(Rp);
%figure(1)
%plot(1:len,Elog,'b-')
%title('Performance Comparison')
%ylabel('Error(dB)')
%xlabel('Number of iterations')
%legend('RTRL')

figure(2)
plot (4000:len,abs(Y_out1(4000:len)'),'r',4000:len,abs(xin(4000:len)),'b')
%plot (1:len,abs(Y_out1(1:len)'),'y')
%figure(2)
%plot(E)
%title('RTRL Algorithm')
%ylabel('Error')
%xlabel('Number of iterations')
%grid;
%legend('CRTRL')

%figure(3)
%plot(1:len,xin,'c+-',1:len,Y_out1,'r.-')
%title('RTRL Algorithm')
%ylabel('Colour')
%xlabel('Number of iterations')
%grid;
%legend('Input Data', 'Predicted Data')

    