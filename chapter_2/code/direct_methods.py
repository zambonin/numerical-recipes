def naive_gaussian_elim(A, b):
    """Gaussian elimination with backward substitution and no pivoting.

    Args:
        A: the matrix of coefficients.
        b: the constant terms.

    Returns:
        x: the solution for the system of equations.
    """
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
    """Gaussian elimination with backward substitution and partial pivoting.

    Args:
        A: the matrix of coefficients.
        b: the constant terms.

    Returns:
        x: the solution for the system of equations.
    """
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


def tridiag_matrix_alg(A, d):
    """
    Check if a matrix is n x n with n > 2 and tridiagonal, and proceed with
    the Thomas algorithm, a simplified form of Gaussian elimination.

    `b`, `a` and `c` are, respectively, the diagonal, subdiagonal and
    superdiagonal of the matrix.

    Args:
        A: the matrix of coefficients.
        d: the constant terms.

    Returns:
        x: the solution for the system of equations.

    Raises:
        ValueError: if the matrix is not tridiagonal, the algorithm shall quit.
    """
    indexes = [-len(A)] + list(range(1, len(A)))
    for row, i in zip(A, indexes):
        if 0 in row[i-1:i+2] or set(row[i+2:] + row[:i-1]) not in ({0}, set()):
            raise ValueError("Not a tridiagonal matrix.")

    b, a, c = [A[-1][-1]], [0], [0]
    for i in range(len(A) - 1):
        b.insert(i, A[i][i])
        a.append(A[i+1][i])
        c.insert(i, A[i][i+1])

    n = len(a)
    for k in range(1, n):
        m = a[k] / b[k-1]
        b[k] -= m * c[k-1]
        d[k] -= m * d[k-1]

    x = [0] * (n-1) + [d[-1] / b[-1]]
    for k in range(n-1, -1, -1):
        x[k-1] = (d[k-1] - c[k-1] * x[k]) / b[k-1]

    return x


def max_residual(x, A, b):
    """Largest residual from solving a system of linear equations.

    Args:
        x: the solution for the system.
        A: the matrix of coefficients.
        b: the constant terms.

    Returns:
        max(r): the maximum value from the list of differences obtained from
                applying the solution to the matrix.
    """
    n = len(A)
    r = [0]*n

    for i in range(n):
        soma = sum(A[i][j] * x[j] for j in range(n))
        r[i] = abs(soma - b[i])

    return max(r)


def partial_pivoting(mat, k, n):
    """Changes rows to achieve an useful matrix setup.

    Args:
        mat: the matrix of coefficients.
        k: the row with an offending value.
        n: the destination row.

    Returns:
        mat: a matrix with switched rows.
    """
    col = [i[k] for i in mat]
    i = col.index(max(col))
    mat[k], mat[i] = mat[i], mat[k]

    return mat
