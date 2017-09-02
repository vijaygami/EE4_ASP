% FUNCTION DRNGD() implements the DRNGD algorithm.
%
% Based on the paper "A class of low complexity and fast convergence algorithms for complex-valued neural networks", % IEEE Workshop on Machine Learning for Signal Processing, pp. 13-22, 2004.
%
% INPUT:
% X_in: input signal which should be scaled according to the dynamic range of nonlinearity 
% desired_x: desired response
% learn_rate: step-size
% L: filter length
% Sample: sample length
% dt_reuse: times of reuse
% AF: AF = 1 sigmoid function; AF = 0 tanh function; 
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
function [e,w,y] = drngd(X_in,desired_x,Sample,learn_rate,L,dt_reuse,W_init,AF);

u = learn_rate;
beta=1;
xin = zeros(L,1) + i*zeros(L,1);
e = zeros(Sample,1) + i*zeros(Sample,1);
w=W_init;


for k=1:Sample
    
    xin = [X_in(k) ; xin(1:L-1)]; 
    
    for t=1:dt_reuse
           
        net(k) = xin.'*w;    

        if AF==1
        sig_function(k) = 1 / ( 1 + exp( -( beta * net(k) ))); % sigmoid function
        sig_function_der(k) = beta * sig_function(k) * ( 1 - sig_function(k) );
        else
        sig_function(k) = tanh(beta * net(k));
        sig_function_der(k) =  beta * ( 1 - (sig_function(k))^2 );
        end
    
        y(k) =sig_function(k);
        e(k) = desired_x(k) - y(k); 
        
        w2 = w + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
        w=w2;
        
        if t==2
            w3 = w2 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w3;
        end
        if t==3
            w4 = w3 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w4;
        end
        if t==4
            w5 = w4 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w5;
        end
        if t==5
            w6 = w5 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w6;
        end
        if t==6
            w7 = w6 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w7;
        end
        if t==7
            w8 = w7 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w8;
        end
         if t==8
            w9 = w8 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w9;
        end
        if t==9
            w10 = w9 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w10;
        end
        if t==10
            w11 = w10 + (u * conj(sig_function_der(k)) * (e(k)) * conj(xin)); 
            w=w11;
        end
       
    end   
    end
