function yPn = fLagrange(n, x, y, xE)
    for k = 1 : length(xE)
        yPn(k) = 0;
        for i = 1 : n+1
            numerador = 1;
            denominador = 1;
            for j = 1 : n+1
                if (j != i)
                    numerador *= (xE(k) - x(j));
                    denominador *= (x(i) - x(j));
                end                        
            end
            produto = numerador / denominador;
            yPn(k) += y(i) * produto;    
        end 
    end
end
