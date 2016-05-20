def max_abs_bound(p):

    def radius_max_abs(p):
        return 1 + (max(p[1:]) / p[0])

    p = [abs(i) for i in p]
    return (radius_max_abs(p), 1 / radius_max_abs(p[::-1]))


def cauchy_bound(p, tol):

    def radius_cauchy(p, tol):
        n = len(p)
        r = abs(p[-1] / p[0]) ** (1/n)
        while True:
            old = r
            r = sum(abs(p[i] / p[0]) * (old ** (n - i))
                    for i in range(1, n)) ** (1/n)
            if (abs(r - old) < tol):
                return r

    return (radius_cauchy(p, tol), 1 / radius_cauchy(p[::-1], tol))


def kojima_bound(p):

    def radius_kojima(p):
        q = [abs(p[i] / p[0]) ** (1/i) for i in range(1, len(p))]
        return sum(sorted(q)[-2:])

    return (radius_kojima(p), 1 / radius_kojima(p[::-1]))


def root_bounds(p, tol):

    bounds = [max_abs_bound(p), cauchy_bound(p, tol), kojima_bound(p)]

    max_radius = min(i[0] for i in bounds)
    min_radius = max(i[1] for i in bounds)
    real = ((max_radius + min_radius) / 2) / (2 ** (1/2))

    return complex(real, -real)
