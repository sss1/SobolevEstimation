function [ f ] = multinormal_dis(mu1,sigma1,mu2,sigma2 )
%output a function which is the square of the difference of two input
%functions.

sigma1_inv = inv(sigma1);
sigma1_det = det(sigma1);
sigma2_inv = inv(sigma2);
sigma2_det = det(sigma2);



    function [y] = pointwise_square(x,y,z)
                
        x = [x;y;z];
        [~,n] = size(x);
        y = zeros(1,n);
        D = 3;
        for i = 1:n
            y(i) = ((2*pi)^(-D/2)*sigma1_det^(-1/2)*exp(-1/2*(x(:,i)-mu1)'*sigma1_inv*(x(:,i)-mu1))-...
                (2*pi)^(-D/2)*sigma2_det^(-1/2)*exp(-1/2*(x(:,i)-mu2)'*sigma2_inv*(x(:,i)-mu2)))^2;
        end
    end

f = @pointwise_square;

end

