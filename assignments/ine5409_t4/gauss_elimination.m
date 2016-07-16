function x = gauss_elimination(n, A)

    for k = 1 : n-1
        A = part_pivoting(n, k, A);
        for i = k+1 : n
            aux = A(i, k) / A(k, k);
            for j = k+1 : n+1
                A(i, j) -= aux * A(k, j);
            end
            A(i, k) = 0;
        end
    end

    if (abs(A(n, n)) < 1e-14)
        if (abs(A(n, n+1)) < 1e-14)
            printf("Sistema indeterminado\n")
            x(n) = 0;
        else
            printf("Sistema impossível\n")
            x(n) = NaN;
        end
    else
        x(n) = A(n, n+1) / A(n, n);
    end

    for i = n-1 : -1 : 1
        soma = 0;
        for j = i+1 : n
            soma += A(i, j) * x(j);
        end
        x(i) = (A(i, n+1) - soma) / A(i, i);
    end

end
