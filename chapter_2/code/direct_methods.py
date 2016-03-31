from numpy import dot
from copy import deepcopy


def naive_gaussian_elim(A, b):
    """Gaussian elimination with backward substitution and no pivoting."""
    A, b, n = deepcopy(A), deepcopy(b), len(A)
    assert len(b) == n, "A and b have different sizes."

    for k in range(n-1):
        for i in range(k+1, n):
            mult = A[i][k] / A[k][k]
            A[i][k] = mult
            for j in range(k+1, n):
                A[i][j] -= mult * A[k][j]
            b[i] -= mult * b[k]

    x = [0]*n
    for i in range(n-1, -1, -1):
        x[i] = (b[i] - dot(A[i][i+1:], x[i+1:])) / A[i][i]

    return x


def max_residual(x, A, b):
    """Largest residual from solving a system of linear equations."""
    n = len(A)
    r = [0]*n

    for i in range(n):
        soma = sum(A[i][j] * x[j] for j in range(n))
        r[i] = abs(soma - b[i])

    return max(r)
