% FUNCTION split_complex_rtrl() implements the Split-Complex RTRL algorithm
%
% Based on the paper of SuLee, "Nonlinear adaptive prediction of complex-valued signals by complex-valued PRNN", IEEE Transactions on Signal Processing, vol 53, no 5, 2005.
% INPUT:
% input: Signal which should be scaled according to the dynamic range of nonlinearity 
%
% OUTPUT:
% Elog: Error in dB
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
function Elog =split_complex_rtrl(input)

%Defining parameters Danilo's book
p=4; %input taps
N=2;%number of neurons
alpha=0.1;%learning rate
len=5000; %length of input sample

%Initialization of weight parameter
w1=zeros(p+N+1,1)+i*zeros(p+N+1,1); % weights of neuron1 +1 for bias yaa...
w2=zeros(p+N+1,1)+i*zeros(p+N+1,1); % weights of neuron2
W=zeros(p+N+1,N)+i*zeros(p+N+1,N); %total weights matrix

%Initialization of weight change
dW1=zeros(1,p+N+1)+i*zeros(1,p+N+1); %weight change +1 for bias also
dW2=zeros(1,p+N+1)+i*zeros(1,p+N+1); %weight change
DW=zeros(N,p+N+1)+i*zeros(N,p+N+1);

%PI for split case 
PII_old=zeros(p+N+1,N,N);%previous sample imaginary
PII_new=zeros(p+N+1,N,N);%current sample imaginary
PRR_old=zeros(p+N+1,N,N);%previous sample real
PRR_new=zeros(p+N+1,N,N);%current sample real
PRI_old=zeros(p+N+1,N,N);%previous sample imaginary
PRI_new=zeros(p+N+1,N,N);%current sample imaginary
PIR_old=zeros(p+N+1,N,N);%previous sample real
PIR_new=zeros(p+N+1,N,N);%current sample real

%Output Matrix
Y_old1=zeros(len+1,1)+i*zeros(len+1,1); % Previous output matrix 1
Y_out1=zeros(len,1)+i*zeros(len,1); % output matrix 1
Y_old2=zeros(len+1,1)+i*zeros(len+1,1); % Previous output matrix 2
Y_out2=zeros(len,1)+i*zeros(len,1); % output matrix 2

%Preparing the basic format fo the input signal
x=zeros(p,1)+i*zeros(p,1);%input signal to the tap
E_dB=zeros(1,len);
Elog=zeros(1,len);

for monte=1:100
    monte;
    
%weight initialization
W_init1real=0.01*rand(p+N+1,1); %initialise the weight1 randomly
W_init1imag=0.01*rand(p+N+1,1); 
w1=W_init1real+i*W_init1imag;
W_init2real=0.01*rand(p+N+1,1); %initialise the weight2 randomly
W_init2imag=0.01*rand(p+N+1,1);
w2=W_init2real+i*W_init2imag;

%Load Data of Complex Colored Input
d=input(1:len);
xin(1)=0;
xin(2:len)=d(1:len-1);

%Activity of the Neurons
for k=1:len
    x=[xin(k);x(1:p-1)]; %input of the system due to input 
    Uin=[Y_old1(k);Y_old2(k);1;x]; %the main input to the system, 1 represents the bias

    %First Neuron Activity
    Vout1=w1.'*Uin;
    Sr1=real(w1).'*real(Uin)-imag(w1).'*imag(Uin);
    Si1=imag(w1).'*real(Uin)+real(w1).'*imag(Uin);
    sig_function1R = 1 ./ ( 1 + exp((-(Sr1))));%Output of the neuron 1 real
    sig_function1I= 1 ./ ( 1 + exp( (-(Si1))));%Output of neuron 1 imag
    sig_function1= sig_function1R + i*sig_function1I;
    sig_function_der1R = sig_function1R .* ( 1 - sig_function1R );%the derivative,f'
    sig_function_der1I = sig_function1I .* ( 1 - sig_function1I );
    Y_out1(k)=sig_function1;%store in the output matrix of 1st neuron
    Y_old1(k+1)=Y_out1(k);
    
    %Output of 1st Neuron
    u1(k)=sig_function1R;
    v1(k)=sig_function1I;
    
    %Second Neuron Activity
    Vout2=w2.'*Uin;
    Sr2=real(w2).'*real(Uin)-imag(w2).'*imag(Uin);
    Si2=imag(w2).'*real(Uin)+real(w2).'*imag(Uin);
    sig_function2R = 1 ./ ( 1 + exp( -(Sr2)));%Output of the neuron 1 real
    sig_function2I= 1 ./ ( 1 + exp( -(Si2)));%Output of neuron 1 imag
    sig_function2= sig_function2R + i*sig_function2I;
    sig_function_der2R = sig_function2R .* ( 1 - sig_function2R );%the derivative,f'
    sig_function_der2I = sig_function2I .* ( 1 - sig_function2I );
    Y_out2(k)=sig_function2;%store in the output matrix of 1st neuron
    Y_old2(k+1)=Y_out2(k);
    
    %Error Calculation
    e(:,k) = d(:,k) - Y_out1(k);
    e_real(k)=real(d(k))-u1(k);
    e_imag(k)=imag(d(k))-v1(k);
    
    %2nd output error calculation
    
    %MSE Error Calculation
    E(k)=(1/2)*((e_real(k)).^2+(e_imag(k)).^2);
    E_dB(k)=10*log10(E(k));%error value at k step
    
    %Matrix ready before calculating Pij
    sig_function_derR=[sig_function_der1R,sig_function_der2R];
    sig_function_derI=[sig_function_der1I,sig_function_der2I];
    %misinterpret formula for splitcomplex; will be back again!
    %Calculating Pij
    for l=1:N%k
        for t=1:N%destination
            for j=1:p+N+1%source
                tempRR=0;
                tempII=0;
                tempRI=0;
                tempIR=0;
                for r=1:N%it is typical to sum both Pij value
                tempRR=tempRR+((real(W(r,l))).*(PRR_old(j,t,r))-imag(W(r,l)).*(PIR_old(j,t,r)));
                tempII=tempII+((real(W(r,l))).*(PII_old(j,t,r))+imag(W(r,l)).*(PRI_old(j,t,r)));
                tempRI=tempRI+((real(W(r,l))).*(PRI_old(j,t,r))-imag(W(r,l)).*(PII_old(j,t,r)));
                tempIR=tempIR+((real(W(r,l))).*(PIR_old(j,t,r))+imag(W(r,l)).*(PRR_old(j,t,r)));
                end
                if l==t
                    tempRR=tempRR+real(Uin(j));
                    tempII=tempII+real(Uin(j));
                    tempRI=tempRI-imag(Uin(j));
                    tempIR=tempIR+imag(Uin(j));
                else
                    tempRR=tempRR;
                    tempII=tempII;
                    tempRI=tempRI;
                    tempIR=tempIR;
                    
                end
                PRR_new(j,t,l)=sig_function_derR(l).*tempRR;
                PII_new(j,t,l)=sig_function_derI(l).*tempII;
                PRI_new(j,t,l)=sig_function_derR(l).*tempRI;
                PIR_new(j,t,l)=sig_function_derI(l).*tempIR;
                
            end
        end
    end
    PRR_old=PRR_new;
    PII_old=PII_new;
    PRI_old=PRI_new;
    PIR_old=PIR_new;
    
    %weight change
    DW=alpha.*(e_real(k).*PRR_new(:,:,1)+e_imag(k).*PIR_new(:,:,1))+i*(e_real(k).*PRI_new(:,:,1)+e_imag(k).*PII_new(:,:,1));
    w1=w1+DW(:,1);
    w2=w2+DW(:,2);
    W=[w1,w2];
 
end%k

Elog=Elog+E_dB;

end%monte
Elog=Elog/monte;


figure(1)
plot(1:len,Elog,'b-')
%title('Split CRTRL Algoritm')
%ylabel('Error(dB)')
%xlabel('Number of iterations')
%grid;
%legend('CRTRL')

%figure(2)
%plot(E)
%title('Split CRTRL Algorithm')
%ylabel('Error')
%xlabel('Number of iterations')
%grid;
%legend('CRTRL')


    
   