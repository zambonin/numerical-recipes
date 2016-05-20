function r = remainders(n, a, xI)
    nI = n;
    for d = 1 : nI
        b(1) = a(1);
        for i = 2 : n+1
            b(i) = a(i) + xI * b(i-1);
        end
        n -= 1;
        a = b;
    end
    r = flip(b);
end
