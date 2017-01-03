function [C erro] = questao1e()

    printf("\n1e)\n");

    xInicial = -1;
    xFinal = 1;
    x = xInicial : 0.01 : xFinal;
    exato = sin(x);

    M = 5;
    C = transpose(coeficientes_exp(M));

    [a b] = fPade(3, 2, C);
    approx = (a(0+1) + a(1+1) .* x + a(2+1) .* (x .^ 2) + a(3+1) .* (x .^ 3)) ./ (b(1) + b(1+1) .* x + b(2+1) .* (x .^ 2));

    erro = max(abs(approx .- exato));

end
