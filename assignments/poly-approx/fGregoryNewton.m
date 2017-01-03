function yPn = fGregoryNewton(n, x, y, xE)
    % Cálculo das diferenças divididas:
    k = 1;
    for i = 1 : n
        DDy(i,k) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    end
    for k = 2 : n
        for i = 1 : n-k+1
            DDy(i,k) = (DDy(i+1, k-1) - DDy(i, k-1)) / (x(i+k) - x(i));
        end
    end
    DDy;

    for indice = 1 : length(xE)  % Para cada ponto de xE
        soma = y(1);
        for k = 1 : n
            produto = 1;
            for j = 1 : k
                produto *= xE(indice) - x(j);
            end
            soma += DDy(1,k) * produto;
        end
        yPn(indice) = soma;
    end
end
