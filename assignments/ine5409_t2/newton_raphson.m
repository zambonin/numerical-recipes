function [x M k difDif] = newton_raphson(n, a, xI, tol, min_val, n_it)
    difAnt = 1;
    difDif = 1;
    k = 0;
    while (difDif > tol && k < n_it)
        k = k + 1;
        r = remainders(n, a, xI);
        M = multiplicity(r, min_val);
        dx = - r(M) / (M * r(M+1));
        x = xI + dx;
        xI = x;
        dif = abs(dx) + abs(r(1));
        difDif = abs(difAnt - dif);
        difAnt = dif;
    end
end
