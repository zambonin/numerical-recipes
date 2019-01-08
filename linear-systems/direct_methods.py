#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""direct_methods.py

A method that allows for computing the solution `x` within a finite number of
operations (in exact arithmetic) and solving the linear system `Ax = b` is
called a direct method.
"""

from __future__ import division


def naive_gaussian_elim(coef, terms):
    """
    Gaussian elimination with backward substitution and no pivoting.

    Args:
        coef:   the matrix of coefficients.
        terms:  the constant terms.

    Returns:
        sol:  the solution for the system of equations.
    """
    assert len(coef) == len(terms), "Matrices have different sizes."
    ext = [[float(k) for k in j] + [terms[i]] for i, j in enumerate(coef)]
    size = len(ext)
    sums = subs = muls = divs = 0

    for i in range(size - 1):
        for j in range(i + 1, size):
            mult = ext[j][i] / ext[i][i]
            divs += 1
            ext[j][i] = mult
            for k in range(i + 1, size):
                ext[j][k] -= mult * ext[i][k]
                subs += 1
                muls += 1
            ext[j][size] -= mult * ext[i][size]
            subs += 1
            muls += 1

    assert ext[size - 1][size - 1] != 0, "No unique solution exists."

    sol = [0] * size
    for i in range(size - 1, -1, -1):
        subs_sum = sum(ext[i][j] * sol[j] for j in range(i + 1, size))
        sums += size - i - 1
        muls += size - i - 1
        sol[i] = (ext[i][size] - subs_sum) / ext[i][i]
        subs += 1
        divs += 1

    print(
        "Contadores separados:\n"
        "Somas: {}    Subtrações: {}\n"
        "Multiplicações: {}    Divisões: {}".format(sums, subs, muls, divs)
    )
    print(
        "Número de operações em ponto flutuante: {}\n".format(
            sums + subs + muls + divs
        )
    )

    return sol
