function A = questao5()

    % 5a)   A interpolação polinomial geralmente é recomendada, pois os pontos
    %       de ancoragem podem ser gerados com precisão.

    printf("\n5{b..d})\n")

    xInicial = 0;
    xFinal = 1;
    n = 2;

    h = (xFinal - xInicial) / n;
    x = xInicial : h : xFinal;
    y = exp(x);
    x = transpose(x);
    y = transpose(y);

    tsis = n+1;
    for i = 1 : tsis
        A(i, tsis) = 1.0;
        for j = tsis-1 : -1 : 1
            A(i, j) = A(i, j+1) * x(i);
        end
        A(i, tsis+1) = y(i);
    end
    C = elimGauss(tsis, A)

    intervalo = xInicial : h*0.05 : xFinal;
    exato = exp(intervalo);

    yPn = fInterpolacaoSistema(n, x, y, intervalo);
    erroPolInt = max(abs(yPn .- exato));

    yPL = fLagrange(n, x, y, intervalo);
    erroLagrange = max(abs(yPL .- exato));

    yPGN = fGregoryNewton(n, x, y, intervalo);
    erroGregNewton = max(abs(yPGN .- exato))

    v = 0.378;
    yV = exp(v);

    yPnv = fInterpolacaoSistema(n, x, y, v);
    erromaxT = max(abs(yPnv .- yV));

    yPLv = fLagrange(n, x, y, v);
    erromaxTL = max(abs(yPLv .- yV));

    yPGNv = fGregoryNewton(n, x, y, v)
    erromaxTGN = max(abs(yPGNv .- yV))

end
