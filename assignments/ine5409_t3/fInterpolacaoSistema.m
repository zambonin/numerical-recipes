function yPn = fInterpolacaoSistema(n, x, y, xE)

    ne = n+1; % número de equações do sistema
    A(ne, ne+1) = 0;
    A(:, ne) = 1;
    A(:, ne+1) = y;
    A(:, ne-1) = x;

    for i = ne-2 : -1 : 1
        A(:, i) = A(:, i+1) .* x;
    end

    yPn = fPn(n, elimGauss(ne, A), xE);

end
