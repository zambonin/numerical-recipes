function [C n] = questao1b()

    printf("\n1b)\n");

    xInicial = -1;
    xFinal = 1;
    tol = sqrt(10) * 1e-2;
    n = 1; % n < 4, termos após o quarto são irrelevantes por conta do erro
    erro = 1;

    while (erro > tol && n < 4)
        n++;
        h = (xFinal - xInicial) / n;
        x = xInicial : h : xFinal;
        serie = [
            0 * (1 / factorial(0));
            1 * (1 / factorial(1));
            0 * (1 / factorial(2));
            -1 * (1 / factorial(3));
            0 * (1 / factorial(5));
            -1 * (1 / factorial(7))
        ];

        for j = 1 : n+1
            yM(j) = 0;
            for i = n : -2 : 0
                yM(j) += serie(i+1) * (x(j))^i;
            end
        end

        exato = sin(x);
        erro = max(abs(yM .- exato));
    end

    C = serie;

end
