function a = poly_adjust(n, m, x, y)

    order = n + 1;
    for i = 1 : order
        for j = 1 : order
            A(i, j) = sum(x(1:m).^(i+j-2));
        end
        b(i) = sum(x(1:m).^(i-1) .* y(1:m));
    end
    a = gauss_elimination(order, [A transpose(b)]);

end
