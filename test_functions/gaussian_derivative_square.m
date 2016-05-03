function [ pdf ] = gaussian_derivative_square(mu,sigma)

    function [y] = gaussian_derivative_square_pdf(x)
        y = (-1/sigma^3/sqrt(2*pi)*(x-mu).*exp(-(x-mu).^2/2/sigma/sigma)).^2;
    end

    pdf = @gaussian_derivative_square_pdf;

end