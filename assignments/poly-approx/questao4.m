function A = questao4()

    % 4a) Experimentando o valor inicial como n = 2, e tomando h como o valor
    % médio de cada intervalo do polinômio de grau n, então

    %     Erro_2 = abs(((-cos(0))*(pi/4)^(2+1)) / (4*(2+1)))
    %            = 0.040373

    % Este erro ainda é maior do que a tolerância sugerida. Portanto, n deve
    % ser incrementado. Para n = 3,

    %     Erro_3 = abs(((sin(pi/2))*(pi/6)^(3+1)) / (4*(3+1)))
    %            = 0.0046976

    % Este erro é válido, pois está abaixo da tolerância sugerida. Então,
    % n = 3 é o grau obtido.

    printf("\n4b)\n")

    xInicial = 0;
    xFinal = pi/2;
    n = 2;

    h = (xFinal - xInicial) / n;
    x = xInicial : h : xFinal;
    y = sin(x);
    x = transpose(x);
    y = transpose(y);

    xE = xInicial : h*0.05 : xFinal;
    yE = sin(xE);

    yPn = fGregoryNewton(n, x, y, xE);
    tsis = n+1;
    for i = 1 : tsis
        A(i, tsis) = 1.0;
        for j = tsis-1 : -1 : 1
            A(i, j) = A(i, j+1) * xE(i);
        end
        A(i, tsis+1) = yPn(i);
    end
    C = elimGauss(tsis, A)

    printf("\n4c)\n")
    erro = max(abs(yPn .- yE))
    erroTaylor = abs(((-cos(0))*(pi/4)^(2+1)) / (4*(2+1)))
    assert (erro < erroTaylor)

    printf("\n4d)\n")
    plot (x, y, "*", xE, yE, "r", xE, yPn, "b")

end
