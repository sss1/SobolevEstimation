function [ pdf ] = uniform_dis(a,b)

    function [y] = unif_pdf(z)
        y = zeros(length(length(z)),1);
        for i = 1:length(z)
            x = z(i);
            if x < a || x > b
                y(i) = 0; 
            else
                y(i) = 1/(b-a);
            end
        end
    end

    pdf = @unif_pdf;

end
