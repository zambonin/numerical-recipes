function xI = root_bounds(n, a)
    % maximum absolute value
        raioMax(1) = 1 + max(abs(a(2:n+1))) / abs(a(1));
        raioMin(1) = 1 / (1 + max(abs(a(1:n))) / abs(a(n+1)));

    % Cauchy bound
        % max:
        rI = abs(a(n+1) / a(1)) ^ (1/n);
        for k = 1 : 30
            soma = abs(a(n+1) / a(1));
            for i = 2 : n
                soma = soma + abs(a(i) / a(1)) * rI ^ (n-i+1);
            end
            raioMax(2) = soma ^ (1/n);
            rI = raioMax(2);
        end

        %min:
        rI = abs(a(1) / a(n+1)) ^ (1/n);
        for k = 1 : 30
            soma = abs(a(1) / a(n+1));
            for i = 2 : n
                soma = soma + abs(a(i) / a(n+1)) * rI ^ (i-1);
            end
            rI = soma ^ (1/n);
        end
        raioMin(2) = 1/rI;

    % Kojima bound
        % max:
        for i = 2 : n+1
            q(i) = abs(a(i) / a(1)) ^ (1 / (i-1));
        end
        ordenado = sort(q);
        max1 = ordenado(n+1);
        max2 = ordenado(n);
        raioMax(3) = max1 + max2;

        % min:
        for i = 2 : n+1
            q(i) = abs(a(n-i+2) / a(n+1)) ^ (1 / (i-1));
        end
        ordenado = sort(q);
        max1 = ordenado(n+1);
        max2 = ordenado(n);
        raioMin(3) = 1 / (max1 + max2);

    raioMedio = ((min(raioMax) + max(raioMin)) / 2);
    a = raioMedio * cos(pi/4);
    xI = complex(a, -a);
end
