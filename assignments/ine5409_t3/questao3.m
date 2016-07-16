function A = questao3()

    np = 100; 
    h = 1 / np; 

    % pontos de ancoragem
    pA = [0, 0];
    pB = [2, 1];
    pC = [8, 0.2];
    pD = [10, 0];

    % pontos de controle
    cA1 = pA + [1, tan(45 * pi/180)];
    cB1 = pB - [.5, tan(0)];
    cB2 = pB + [1, tan(-15 * pi/180)];
    cC1 = pC - [1, tan(-8 * pi/180)];
    cC2 = pC + [1, tan(-8 * pi/180)];
    cD1 = pD - [1, tan(-5 * pi/180)];

    bezierAB = bezier(np, pA, cA1, pB, cB1);
    bezierBC = bezier(np, pB, cB2, pC, cC1);
    bezierCD = bezier(np, pC, cC2, pD, cD1);

    points = [pA; pB; pC; pD];
    control = [cA1; cB1; cB2; cC1; cC2; cD1];
    final = [bezierAB; bezierBC; bezierCD];

    plot(
        final(:, 1), final(:, 2), 'r',
        points(:, 1), points(:, 2), '*b',
        control(:, 1), control(:, 2), '*k'
    )

end

function x = bezier(np, p1, c1, p2, c2)

    h = 1 / np;

    c = 3 * (c1 - p1);
    b = 3 * (c2 - c1) - c;
    a = (p2 - p1) - (c + b);

    t = 0;
    for i = 1 : np+1
        x(i, :) = p1 + t * (c + t * (b + t * a));
        t += h;
    end

end
