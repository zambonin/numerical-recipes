function [x M iters] = orig_roots(a)
    tol = 1e-12;
    n_it = 120;
    min_val = 1e-4;
    n = length(a) - 1;
    nI = n;
    aI = a;
    i_raiz = 0;
    it_1 = 0;
    it_2 = 0;
    iters = 0;

    while n > 0
        i_raiz += 1;
        printf("Raiz %d\n", i_raiz)

        printf("    a)  Valor inicial: "), xI = root_bounds(n, a)

        [x(i_raiz) M(i_raiz) k difDif] = newton_raphson(n, a, xI, tol, min_val,
                                                        n_it);
        printf("    c)  Iterações acumuladas na primeira chamada do Newton-Raphson: ")
        it_1 += k

        [x(i_raiz) M(i_raiz) k difDif] = newton_raphson(nI, aI, x(i_raiz), tol,
                                                        min_val, n_it);
        printf("    d)  Iterações acumuladas na segunda chamada do Newton-Raphson: ")
        it_2 += k

        printf("    e)  Critério de parada atingido: "), difDif
        [n, a] = poly_div(n, a, x(i_raiz), M(i_raiz));
    end
    iters = (it_1 + it_2);
end
