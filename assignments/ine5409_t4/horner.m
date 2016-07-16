function y = horner(n, a, xi)

    for k = 1 : length(xi)
        y(k) = a(n+1);
        for i = n:-1:1
            y(k) = a(i) + y(k) * xi(k);
        end
    end

end
