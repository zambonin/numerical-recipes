from copy import deepcopy
from sys import float_info


def jacobi_method(A, b, n_it, tol=float_info.epsilon):
    """Jacobi iterative method used to solve sparse matrices.

    Args:
        A: the matrix of coefficients.
        b: the constant terms.
        n_it: number of iterations.
        tol: metric to check if the precision of the values are good enough
            and prevent further iterations. Default: double precision
            machine epsilon (2^-52, approx. 2.22e-16).

    Returns:
        x: the solution for the system of equations.
    """
    assert len(A) == len(b), "A and b have different sizes."
    ext = [j + [b[i]] for i, j in enumerate(A)]
    n = len(ext)

    x = [0]*n
    for k in range(n_it):
        prev = deepcopy(x)
        for i in range(n):
            summation = sum(ext[i][j] * prev[j] for j in range(n) if j != i)
            x[i] = (ext[i][n] - summation) / ext[i][i]
        if abs(max(prev[i] - x[i] for i in range(len(x)))) < tol:
            print("Maximum tolerance exceeded at iteration {}.".format(k+1))
            break

    return x


def SOR(A, b, n_it, relax, tol=float_info.epsilon):
    """Successive over-relaxation is a variant of Gauss-Seidel, where a
    chosen relaxation factor may result in faster or slower convergence.

    Args:
        A: the matrix of coefficients.
        b: the constant terms.
        n_it: number of iterations.q
        relax: relaxation factor where 0 < ω < 2. Conventional Gauss-Seidel is
            achieved with ω = 1.
        tol: metric to check if the precision of the values are good enough
            and prevent further iterations. Default: double precision
            machine epsilon (2^-52, approx. 2.22e-16).

    Returns:
        x: the solution for the system of equations.
    """
    assert len(A) == len(b), "A and b have different sizes."
    ext = [j + [b[i]] for i, j in enumerate(A)]
    n = len(ext)

    x = [0]*n
    for k in range(n_it):
        prev = deepcopy(x)
        for i in range(n):
            summation = sum(ext[i][j] * x[j] for j in range(n) if j != i)
            x[i] += relax * (((ext[i][n] - summation) / ext[i][i]) - x[i])
        if abs(max(prev[i] - x[i] for i in range(len(x)))) < tol:
            print("Maximum tolerance exceeded at iteration {}.".format(k+1))
            break

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
        b[k] -= m*c[k-1]
        d[k] -= m*d[k-1]

    x = [0] * (n-1) + [d[-1] / b[-1]]
    for k in range(n-1, -1, -1):
        x[k-1] = (d[k-1] - c[k-1] * x[k]) / b[k-1]

    return x
