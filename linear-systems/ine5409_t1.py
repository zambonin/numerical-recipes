#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ine5409_t1.py

Code for the first assignment, which solves sparse linear systems with
direct and iterative numerical methods.
"""

from direct_methods import naive_gaussian_elim
from iterative_methods import (successive_over_relaxation,
                               check_diagonal_dominance)


def create_matrix():
    """Creates a system of 100 linear equations for the assignment."""
    size = 100
    coef, terms = [[0] * size for i in range(size)], [0] * size

    coef[0][0], coef[0][1], terms[0] = 1, 1, 1.5

    for i in range(1, size // 2):
        coef[i][i - 1], coef[i][i], coef[i][i + 1] = 1, 3, 1
        coef[i][i + (size // 2)], terms[i] = 1, 1

    for i in range(size // 2, size - 1):
        coef[i][i - size // 2], coef[i][i - 1], coef[i][i] = 1, 1, 3
        coef[i][i + 1], terms[i] = 1, 2

    coef[-1][-2], coef[-1][-1], terms[-1] = 1, 1, 3

    return coef, terms


def main():
    """Executes the numerical methods and prints solutions."""
    coef, terms = create_matrix()

    # questão 1
    q_1 = naive_gaussian_elim(coef, terms)
    print("Solução por eliminação Gaussiana: {}\n".format(q_1))

    # questão 2
    print("Dominância diagonal: {}\n".format(check_diagonal_dominance(coef)))

    q_2 = successive_over_relaxation(coef, terms, 500, 1.879, 1e-4)
    print("Solução por SOR: {}\n".format(q_2))

    trunc = max(abs((i - j) / j) for i, j in zip(q_2, q_1))
    print("Maior erro de truncamento: {}".format(trunc))


if __name__ == "__main__":
    main()
