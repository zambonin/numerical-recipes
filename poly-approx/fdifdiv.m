function difdiv1 = fdifdiv(n, x, y)
    k = 1;

    for i = 1 : n
        difdiv(i, 1) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    end

    for k = 2 : n
        for i = 1 : n+1-k
            difdiv(i, k) = (difdiv(i+1, k-1) - difdiv(i, k-1)) / (x(i+k) - x(i));
        end
    end

    difdiv1 = difdiv(1, :);
end
