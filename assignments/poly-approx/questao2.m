function A = questao2()

    x = [0, 1, 2, 3];
    y = [-3, -2, 4, 0];
    n = length(x)-1;

    printf("\n2a)\n");

    for i = 1 : n
        h(i) = x(i+1) - x(i);
    end

    [a, b, c, d, S] = fSplinea(n, x, y, h);
    a = a(2)
    b = b(2)
    c = c(2)
    d = d(2)
    S

    printf("2b)\n")

    [a, b, c, d, S] = fSplineb(n, x, y, h);
    a = transpose(a)
    b = transpose(b)
    c = transpose(c)
    d = transpose(d)
    S = transpose(S)

    printf("2c)\n")

    xp = 1.3;
    is = 0;
    for i = 1 : n
        if xp >= x(i)
            is = i;
        end
    end
    is

    printf("\n2d)\n")
    l = xp - x(2);
    yp = a(2)*(l^3) + b(2)*(l^2) + c(2)*(l) + d(2);
    yp

    % 2e)
    np = 10;
    xpp = [];
    ypp = [];
    for i = 1 : n
        xp = x(i) : (x(i+1)-x(i))/np : x(i+1);
            for k = 1 : np+1
                l = xp(k) - x(i);
               yp(k) = a(i)*(l^3) + b(i)*(l^2) + c(i)*(l) + d(i);
            end
            xpp=[xpp xp];
            ypp=[ypp yp];
    end

    plot(x, y, '*' , xpp, ypp, 'r', 5, 0)

end