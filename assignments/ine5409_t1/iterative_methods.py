# -*- coding: utf-8 -*-

from copy import deepcopy
from sys import float_info


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
    n = len(A)
    sums = subs = muls = divs = 0

    x = [0]*n
    for k in range(n_it):
        prev = deepcopy(x)

        i = 0
        x[i] = prev[i] + relax * (((b[i] - x[i+1]) / A[i][i]) - prev[i])

        for i in range(1, n//2):
            sumx50 = x[i-1] + x[i+1] + x[i + (n//2)]
            x[i] = prev[i] + relax * (((b[i] - sumx50) / A[i][i]) - prev[i])

        for i in range(n//2, n-1):
            sumx99 = x[i - (n//2)] + x[i-1] + x[i+1]
            x[i] = prev[i] + relax * (((b[i] - sumx99) / A[i][i]) - prev[i])

        i = n - 1
        x[-1] = prev[-1] + relax * (((b[-1] - x[i-1]) / A[-1][-1]) - prev[-1])

        sums += 2 + 3 * len(range(n-2))
        muls += len(range(n))
        subs += 2 * len(range(n))
        divs += len(range(n))

        if max(abs((i - j) / j) for i, j in zip(prev, x)) < tol:
            print("Maximum tolerance exceeded at iteration {}.".format(k+1))
            break

    print("Contadores separados:\n"
          "Somas: {}    Subtrações: {}\n"
          "Multiplicações: {}    Divisões: {}"
          .format(sums, subs, muls, divs))
    print("Número de operações em ponto flutuante: {}\n".format(
        sums + subs + muls + divs))

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
