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


def SOR(A, b, n_it, relax=1, tol=float_info.epsilon):
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


def check_diagonal_dominance(mat):
    """Check if a matrix is diagonally dominant.

    Args:
        mat: the matrix of coefficients.

    Returns:
        The all operator returns True if all elements of a list are
        also True. Hence, if any element fails to be dominant over
        its row, the matrix itself cannot be diagonally dominant.
    """
    return all(abs(mat[i][i]) > sum(abs(k) for k in j) - abs(mat[i][i])
               for i, j in enumerate(mat))
