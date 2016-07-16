function yPn = fPn(n, a, xE)
    % Briot-Ruffini
    for k = 1 : length(xE)
        b(1) = a(1);
        for i = 2 : n+1
            b(i) = a(i) + xE(k) * b(i-1);
        end
        yPn(k) = b(n+1);
    end
end
