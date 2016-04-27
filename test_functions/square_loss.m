function [ f ] = square_loss(f_1,f_2 )
%output a function which is the square of the difference of two input
%functions.

f = @(x) (f_1(x)-f_2(x)).^2;


end

