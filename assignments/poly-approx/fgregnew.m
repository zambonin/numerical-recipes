function f = fgregnew(xp, n, x, y, difdiv1)
    f = y(1);
    for k = 1 : n
        prod = 1;
        for j = 1 : k
            prod = prod * (xp - x(j));
        end
        f = f + difdiv1(k) * prod;
    end
end
