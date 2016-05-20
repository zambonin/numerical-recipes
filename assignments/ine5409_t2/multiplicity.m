function M = multiplicity(r, min_val)
    rLimite = min_val * r(length(r));
    M = 1;
    somaRestos = abs(r(1)) + abs(r(2));
    while somaRestos < rLimite
        M = M + 1;
        somaRestos = somaRestos + abs(r(M+1));
    end
end
