function erromax = questao6b()

    printf("\n6b)\n");

    xInicial = 0;
    xFinal = pi/2;
    tol = sqrt(10) * 1e-6;
    n_it = 1;
    n = 1;
    erromax = 1;

    while (erromax > tol && n_it < 30)
        n_it++;
        n++;
        h = (xFinal - xInicial) / n;
        x = xInicial : h : xFinal;
        y = cos(x);

        difdiv1 = fdifdiv(n, x, y);
        np = 20*n;
        hp = (x(n+1) - x(1)) / np;
        xp = x(1) : hp : x(n+1);
        ye = cos(xp);

        for i = 1 : np+1
            yip(i) = fgregnew(xp(i), n, x, y, difdiv1);
        end

        erro = abs(yip .- ye);
        erromax = max(erro);
        plot(xp, erro, "r;Erro entre Pn(x) e f(x);")
    end

end
