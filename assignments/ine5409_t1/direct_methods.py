# -*- coding: utf-8 -*-


def naive_gaussian_elim(A, b):
    """Gaussian elimination with backward substitution and no pivoting.
    Args:
        A: the matrix of coefficients.
        b: the constant terms.
    Returns:
        x: the solution for the system of equations.
    """
    assert len(A) == len(b), "A and b have different sizes."
    ext = [[float(k) for k in j] + [b[i]] for i, j in enumerate(A)]
    n = len(ext)
    sums = subs = muls = divs = 0

    for i in range(n-1):
        for j in range(i+1, n):
            mult = ext[j][i] / ext[i][i]
            divs += 1
            ext[j][i] = mult
            for k in range(i+1, n):
                ext[j][k] -= mult * ext[i][k]
                subs += 1
                muls += 1
            ext[j][n] -= mult * ext[i][n]
            subs += 1
            muls += 1

    assert ext[n-1][n-1] != 0, "No unique solution exists."

    x = [0]*n
    for i in range(n-1, -1, -1):
        subs_sum = sum(ext[i][j] * x[j] for j in range(i+1, n))
        sums += len(range(i+1, n))
        muls += len(range(i+1, n))
        x[i] = (ext[i][n] - subs_sum) / ext[i][i]
        subs += 1
        divs += 1

    print("Contadores separados:\n"
          "Somas: {}    Subtrações: {}\n"
          "Multiplicações: {}    Divisões: {}"
          .format(sums, subs, muls, divs))
    print("Número de operações em ponto flutuante: {}\n".format(
        sums + subs + muls + divs))

    return x
