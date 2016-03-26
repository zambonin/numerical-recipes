def quad_eq(a, b, c):
    d = (b**2 - 4*a*c)**(1/2)
    return (- b / (2*a)) if d == 0 else ((- b + d) / (2*a), (- b - d) / (2*a))


def rat_quad_eq(a, b, c):
    d = (b**2 - 4*a*c)**(1/2)
    return (- b / (2*a)) if d == 0 else ((- 2*c) / (b + d), (- 2*c) / (b - d))


print("x1 = {:.25f} x2 = {:.25f}".format(*quad_eq(1, 62.1, 1)))
print("x1 = {:.25f} x2 = {:.25f}".format(*rat_quad_eq(1, 62.1, 1)))
print("x1 = {} {: <20} x2 = {}".format(round(-0.0161072, 4), "",
                                       round(-62.0839, 2)))
print("x1 = {} {: <17} x2 = {}".format(-0.0161072, "", -62.0839))
