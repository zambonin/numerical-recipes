function Tn = trapezoidal_rule(n, a, b)

    h = (b - a) / n;
    x = a : h : b;
    y = exp(-x .^ 2);

    Tn = h/2 * (y(1) + 2 * sum(y(2:n)) + y(n+1));

end
