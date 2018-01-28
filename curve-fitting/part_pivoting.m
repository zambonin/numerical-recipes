function A = part_pivoting(n, k, A)

    [MAX, i] = max(abs(A(k:n, k)));
    i += k - 1;
    aux = A(k, :);
    A(k, :) = A(i, :);
    A(i, :) = aux;

end
