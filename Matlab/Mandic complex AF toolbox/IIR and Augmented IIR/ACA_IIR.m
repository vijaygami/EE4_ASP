% FUNCTION ACA_IIR() provides an estimate of a filter's coefficients using a stochastic gradient method.
%
% Based on the paper of C. Cheong Took and D. Mandic, "Adaptive IIR Filtering of Noncircular Complex Signals", IEEE Transactions on Signal Processing, to be published, 2009.
%
% INPUT:
% x: filter input [(N+1) x 1]
% d: desired response [(N+1) x 1]
% M: number of taps
% mu: step-size
%
% OUTPUT:
% a, b, h, g : estimated taps [M x 1]
% e: estimation error for each tap [1 x n]
% y: filter output
% y(n)= dot(a,y) + dot(b,x) + dot(g,y^*) + dot(h,x^*) ;
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
function [a,b,h,g,e,y] = ACA_IIR(x,d,M,mu)

L=length(x); % length of data
N=L-1 ; 
StepA=1; % one-step ahead prediction

for stepAhead=1:StepA

    y = (zeros(1,L) + i*zeros(1,L))';
    e=y;

    a = transpose(zeros(1,M-1) + i*zeros(1,M-1));
    b = transpose(zeros(1,M) + i*zeros(1,M));
    h = transpose(zeros(1,M) + i*zeros(1,M));
    g = transpose(zeros(1,M-1) + i*zeros(1,M-1));

    psia= transpose(zeros(1,M-1) + i*zeros(1,M-1));
    psib = transpose(zeros(1,M) + i*zeros(1,M));
    psig= transpose(zeros(1,M-1) + i*zeros(1,M-1));
    psih = transpose(zeros(1,M) + i*zeros(1,M));

    phia= transpose(zeros(1,M-1) + i*zeros(1,M-1));
    phib = transpose(zeros(1,M) + i*zeros(1,M));
    phig= transpose(zeros(1,M-1) + i*zeros(1,M-1));
    phih = transpose(zeros(1,M) + i*zeros(1,M));


    %begin of algorithm
    for n = M : N
        if (n+1+stepAhead)>N
            break
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INPUTS
        u = x(n:-1:n-M+1) ;
        conju=conj(u);

        if n==M
            y(1:M) = u(1:M);
        end

        yy = y(n-1:-1:n-M+1) ;
        conjyy = conj(yy);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        y(n+stepAhead)= dot(a,yy)     + dot(b,u) +...
            dot(g,conjyy) + dot(h,conju) ;


        e(n) = d(n+stepAhead) - y(n+stepAhead) ;
        conje = conj(e(n));

        %%%%%%%%%%%%Updates of sensitivities or gradients%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n==M
            phia = conjyy;
            phib = conju;
            phig = yy;
            phih = u;

        else
            phia = [conjyy(1) + dot(conj(a),phia)           +  dot(conj(g),psia);          phia(1:end-1)] ;
            phib = [conju(1)  + dot(conj(a),phib(1:end-1))  +  dot(conj(g),psib(1:end-1)); phib(1:end-1)];
            phig = [yy(1)     + dot(conj(a),phig)           +  dot(conj(g),psig);          phig(1:end-1)];
            phih = [u(1)      + dot(conj(a),phih(1:end-1))  +  dot(conj(g),psih(1:end-1)); phih(1:end-1)];

            psia= [dot(a,psia)          + dot(g,phia);          psia(1:end-1)];
            psib= [dot(a,psib(1:end-1)) + dot(g,phib(1:end-1)); psib(1:end-1)];
            psig= [dot(a,psig)          + dot(g,phig);          psig(1:end-1)];
            psih= [dot(a,psih(1:end-1)) + dot(g,phih(1:end-1)); psih(1:end-1)];
        end


        %%%%%%%%%%%%Updates of coeffs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n< 101 % only 100 first samples to train, then use same coefficients
            mua = mu;
            mub = 5*mu; 
            muh = (10^-1)*mu; 
            mug = mu;
        else
            mua=0;
            mub=0;
            muh = 0;
            mug = 0;
        end

        a = a + mua*( e(n)*phia + psia*conje );
        b = b + mub*( e(n)*phib + psib*conje);
        h = h + muh*( e(n)*phih + psih*conje );
        g = g + mug*( e(n)*phig + psig*conje);
    end

end


%Plotting of graphs
figure,subplot(211), plot(real(y),'k-.')
hold on,plot(real(x),'r')
legend('ACA IIR','Actual'), grid
axis([1 length(x) min(real(x)) max(real(x))]), xlabel('Time (samples)')
subplot(212), plot(imag(y),'k-.') 
hold on,plot(imag(x),'r'), grid
axis([1 length(x) min(imag(x)) max(imag(x))]), xlabel('Time (samples)')


