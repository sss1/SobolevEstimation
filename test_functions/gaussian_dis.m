function [ pdf ] = gaussian_dis(mu,sigma)
%center at middle point of a and b
D = size(mu,1);

sigma_inv = inv(sigma);
sigma_det = det(sigma);

    function [y] = gaussian_pdf(x)
        [~,n] = size(x);
        y = zeros(1,n);
        for i = 1:n
            y(i) = (2*pi)^(-D/2)*sigma_det^(-1/2)*exp(-1/2*(x(:,i)-mu)'*sigma_inv*(x(:,i)-mu));
        end
    end

    pdf = @gaussian_pdf;

end