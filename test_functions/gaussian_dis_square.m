function [ pdf ] = gaussian_dis_square(mu,sigma)

    function [y] = gaussian_square_pdf(x)
        y = (1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/sigma/sigma/2)).^2;
    end

    pdf = @gaussian_square_pdf;

end