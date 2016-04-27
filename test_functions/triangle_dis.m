function [ pdf ] = triangle_dis(a,b)
%center at middle point of a and b

    function [y] = triangle_pdf(z)
        y = zeros(length(length(z)),1);
        for i = 1:length(z);
            x = z(i);
            if x < a || x > b
                y(i) = 0; 
            elseif x < (a+b) / 2
                y(i) = (x - a).*(4/(b-a));
            else
                y(i) = 2/(b-a) - (x - (a+b)/2).*(4/(b-a));
            end
        end
        
    end

    pdf = @triangle_pdf;

end

