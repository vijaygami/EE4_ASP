function [ xp ] = xpast( x, n, p )
% Function taken from my (Vijay Gami) ASP courcework year 3.
% returns the current and past p outputs of x as a column vector 
%  [x(n), x[n-1], ..., x(n-p)]^T (Transposed)

    for j=1:p+1,
        if((n-j) >= 0) xp(j,1) = x(n-j+1);    %p past inputs
        else xp(j,1) =0;                      %x is casual, so x[n]=0 for n<0.  
        end
    end

end

