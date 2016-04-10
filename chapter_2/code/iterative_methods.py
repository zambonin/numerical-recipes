from copy import deepcopy
from sys import float_info


def jacobi_method(A, b, n_it, tol=float_info.epsilon):
    """Jacobi iterative method used to solve sparse matrices."""
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


def gauss_seidel(A, b, n_it, tol=float_info.epsilon):
    """Gauss-Seidel iterative technique, a possible improvement for Jacobi."""
    assert len(A) == len(b), "A and b have different sizes."
    ext = [j + [b[i]] for i, j in enumerate(A)]
    n = len(ext)

    x = [0]*n
    for k in range(n_it):
        prev = deepcopy(x)
        for i in range(n):
            summation = sum(ext[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (ext[i][n] - summation) / ext[i][i]
        if abs(max(prev[i] - x[i] for i in range(len(x)))) < tol:
            print("Maximum tolerance exceeded at iteration {}.".format(k+1))
            break

    return x


def SOR(A, b, n_it, relax, tol=float_info.epsilon):
    """Successive over-relaxation is a variant of Gauss-Seidel."""
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


A = [[10, -1, 2, 0], [-1, 11, -1, 3], [2, -1, 10, -1], [0, 3, -1, 8]]
b = [6, 25, -11, 15]
print(jacobi_method(A, b, 50))
print(gauss_seidel(A, b, 50))
print(SOR(A, b, 50, 0.75))
