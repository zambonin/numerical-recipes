#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""iterative_methods.py

An iterative method is an algorithm that generates a sequence of improving
approximate solutions for a class of problems (in this case, solving linear
systems), in which the `n`-th approximation is derived from the previous ones.
"""

from copy import deepcopy
from sys import float_info


def successive_over_relaxation(coef, terms, n_it, relax=1,
                               tol=float_info.epsilon):
    """
    Successive over-relaxation is a variant of Gauss-Seidel, where a chosen
    relaxation factor may result in faster or slower convergence.

    Args:
        coef:   the matrix of coefficients.
        terms:  the constant terms.
        n_it:   number of iterations.
        relax:  relaxation factor where 0 < ω < 2. Conventional Gauss-Seidel
                is achieved with ω = 1.
        tol:    metric to check if the precision of the values are good enough
                and prevent further iterations. Default: double precision
                machine epsilon (2^-52 or approximately 2.22e-16).

    Returns:
        sol:    the solution for the system of equations.
    """
    assert len(coef) == len(terms), "Matrices have different sizes."
    size = len(coef)
    sums = subs = muls = divs = 0

    sol = [0] * size
    for k in range(n_it):
        prev = deepcopy(sol)

        i = 0
        sol[i] = prev[i] \
                + relax * (((terms[i] - sol[i + 1]) / coef[i][i]) - prev[i])

        for i in range(1, size // 2):
            partial = sol[i - 1] + sol[i + 1] + sol[i + (size // 2)]
            sol[i] = prev[i] \
                    + relax * (((terms[i] - partial) / coef[i][i]) - prev[i])

        for i in range(size // 2, size - 1):
            partial = sol[i - (size // 2)] + sol[i - 1] + sol[i + 1]
            sol[i] = prev[i] \
                    + relax * (((terms[i] - partial) / coef[i][i]) - prev[i])

        i = size - 1
        sol[i] = prev[i] \
                + relax * (((terms[i] - sol[i - 1]) / coef[i][i]) - prev[i])

        sums += 2 + 3 * len(range(size - 2))
        muls += len(range(size))
        subs += 2 * len(range(size))
        divs += len(range(size))

        if max(abs((i - j) / j) for i, j in zip(prev, sol)) < tol:
            print("Maximum tolerance exceeded at iteration {}.".format(k + 1))
            break

    print("Contadores separados:\n"
          "Somas: {}    Subtrações: {}\n"
          "Multiplicações: {}    Divisões: {}".format(sums, subs, muls, divs))
    print("Número de operações em ponto flutuante: {}\n".format(
        sums + subs + muls + divs))

    return sol


def check_diagonal_dominance(mat):
    """Check if a matrix is diagonally dominant.

    Args:
        mat:    the matrix of coefficients.

    Returns:
        The all operator returns True if all elements of a list are
        also True. Hence, if any element fails to be dominant over
        its row, the matrix itself cannot be diagonally dominant.
    """
    return all(abs(mat[i][i]) > sum(abs(k) for k in j) - abs(mat[i][i])
               for i, j in enumerate(mat))
