function [ pdf ] = gaussian_dis(mu,sigma)
%center at middle point of a and b

    function [y] = gaussian_pdf(x)
        y = 1/sqrt(2*pi)/sigma*exp(-(x-mu).^2/sigma/sigma/2);
    end

    pdf = @gaussian_pdf;

end