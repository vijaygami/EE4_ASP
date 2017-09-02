function [ c_coef ] = circ_coef( x )
% Returns the circularity coefficinet of the input data

 cv=cov(x, conj(x)); cv=cv(1,2);
 c_coef = cv/(sqrt(var(x)).*sqrt(var(conj(x))));
 
 %c_coef = mean(x.*x)./mean(x.*conj(x));

end

