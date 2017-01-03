function Sn = simpson_rule(n, a, b)

    h = (b - a) / n;
    x = a : h : b;
    y = exp(-x .^ 2);

    Sn = h/3 * (y(1) + 4 * sum(y(2:2:n)) + 2 * sum(y(3:2:n-1)) + y(n+1));

end
