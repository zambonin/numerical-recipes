function [C n] = questao1a()

    printf("1a)\n");

    xInicial = -1;
    xFinal = 1;
    tol = sqrt(10) * 1e-2;
    n_it = 1;
    n = 1; % grau de P_n(x)
    erro = 1;

    while (erro > tol && n_it < 30)
        n_it++;
        n++;
        h = (xFinal - xInicial) / n;
        x = xInicial : h : xFinal;
        y = sin(x);

        tsis = n+1;
        for i = 1 : tsis
            A(i, tsis) = 1.0;
            for j = tsis-1 : -1 : 1
                A(i, j) = A(i, j+1) * x(i);
            end
            A(i, tsis+1) = y(i);
        end

        C = elimGauss(tsis, A); % coeficientes do polinômio interpolador
        intervalo = xInicial : h*0.05 : xFinal;

        x = transpose(x);
        y = transpose(y);

        exato = sin(intervalo);
        approx = fGregoryNewton(n, x, y, intervalo);
        erro = max(abs(approx .- exato));
    end

end
