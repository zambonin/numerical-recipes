# -*- coding: utf-8 -*-

#####
# Relatório:
# https://github.com/zambonin/UFSC-INE5409/blob/master/assignments/ine5409_t1/ine5409_t1.pdf
#####

from direct_methods import naive_gaussian_elim
from iterative_methods import SOR, check_diagonal_dominance


def create_matrix():

    n = 100
    A = [[0]*100 for i in range(n)]
    b = [0]*n

    A[0][0], A[0][1], b[0] = 1, 1, 1.5

    for i in range(1, n//2):
        A[i][i-1], A[i][i], A[i][i+1] = 1, 3, 1
        A[i][i + (n//2)], b[i] = 1, 1

    for i in range(n//2, n-1):
        A[i][i - n//2], A[i][i-1], A[i][i] = 1, 1, 3
        A[i][i+1], b[i] = 1, 2

    A[-1][-2], A[-1][-1], b[-1] = 1, 1, 3

    return A, b


if __name__ == "__main__":
    A, b = create_matrix()

    # questão 1
    w = naive_gaussian_elim(A, b)
    print("Solução por eliminação Gaussiana: {}\n".format(w))

    # questão 2
    print("Dominância diagonal: {}\n".format(check_diagonal_dominance(A)))

    # r = [i/1000 for i in range(1800, 1900)]
    # for i in r:
    #     x = SOR(A, b, 1000, i)
    #     print("{:.3f}".format(i), x[0] + x[1])

    z = SOR(A, b, 500, 1.879, 1e-4)
    print("Solução por SOR: {}\n".format(z))

    trunc = max(abs((i - j) / j) for i, j in zip(z, w))
    print("Maior erro de truncamento: {}".format(trunc))
