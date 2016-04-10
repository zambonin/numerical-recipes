def naive_gaussian_elim(A, b):
    """Gaussian elimination with backward substitution and no pivoting."""
    assert len(A) == len(b), "A and b have different sizes."
    ext = [j + [b[i]] for i, j in enumerate(A)]
    n = len(ext)

    for i in range(n-1):
        for j in range(i+1, n):
            mult = ext[j][i] / ext[i][i]
            ext[j][i] = mult
            for k in range(i+1, n):
                ext[j][k] -= mult * ext[i][k]
            ext[j][n] -= mult * ext[i][n]

    assert ext[n-1][n-1] != 0, "No unique solution exists."

    x = [0]*n
    for i in range(n-1, -1, -1):
        subs_sum = sum(ext[i][j] * x[j] for j in range(i+1, n))
        x[i] = (ext[i][n] - subs_sum) / ext[i][i]

    return x


def part_pivot_gaussian_elim(A, b):
    """Gaussian elimination with backward substitution and partial pivoting."""
    assert len(A) == len(b), "A and b have different sizes."
    ext = [j + [b[i]] for i, j in enumerate(A)]
    n = len(ext)

    for i in range(n-1):
        ext = partial_pivoting(ext, i, n)
        for j in range(i+1, n):
            mult = ext[j][i] / ext[i][i]
            ext[j][i] = mult
            for k in range(i+1, n):
                ext[j][k] -= mult * ext[i][k]
            ext[j][n] -= mult * ext[i][n]

    assert ext[n-1][n-1] != 0, "No unique solution exists."

    x = [0]*n
    for i in range(n-1, -1, -1):
        subs_sum = sum(ext[i][j] * x[j] for j in range(i+1, n))
        x[i] = (ext[i][n] - subs_sum) / ext[i][i]

    return x


def max_residual(x, A, b):
    """Largest residual from solving a system of linear equations."""
    n = len(A)
    r = [0]*n

    for i in range(n):
        soma = sum(A[i][j] * x[j] for j in range(n))
        r[i] = abs(soma - b[i])

    return max(r)


def partial_pivoting(mat, k, n):
    """Changes rows to achieve an useful matrix setup."""
    col = [i[k] for i in mat]
    i = col.index(max(col))
    mat[k], mat[i] = mat[i], mat[k]

    return mat