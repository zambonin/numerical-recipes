function [n a] = poly_div(n, a, x, M)
    for d = 1 : M
        b(1) = a(1);
        for i = 2 : n+1
            b(i) = a(i) + x * b(i-1);
        end
        n -= 1;
        a = b;
    end
    aux = a(1:n+1);
    a = aux;
end
