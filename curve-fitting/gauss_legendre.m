function [gm, c, x, y] = gauss_legendre(m, a, b)

    [C, t] = gl_coef();

    temp = 0;
    for k = 1 : m
        x(k) = 1/2 * ((b - a) * t(m, k) + (b + a));
        y(k) = exp(-x(k) .^ 2);
        temp += C(m,k) * y(k);
    end

    gm = 1/2 * (b - a) * temp;
    y *= 2 / sqrt(pi);

    % fits a polynomial of degree m - 1 on points (x, y)
    c = polyfit(x, y, m - 1);

end
