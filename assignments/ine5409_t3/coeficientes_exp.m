function c = coeficientes_exp(n)

	for i = 1 : n
        j = 2*i;
        c(j-1) = 0;
        c(j) = (-1)^(i-1) / factorial(2*i-1);
	end

end
