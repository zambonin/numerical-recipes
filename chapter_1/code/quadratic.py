# Numerical Analysis (Burden & Faires, 9th ed., p. 25)
# Demonstration of how subtracting floating point numbers affects greatly the
# final result. The methods purposely round numbers to show errors that may
# appear, and thus, shall not be used for real solving purposes.


def quad_eq(a, b, c):
    """Quadratic formula simple implementation."""
    d = round((b**2 - 4 * a * c)**(1.0 / 2), 2)
    return -b / (2 * a) if d == 0 else ((-b + d) / (2 * a), (-b - d) / (2 * a))


def rat_quad_eq(a, b, c):
    """Rationalized quadratic formula. It should provide
    further accuracy with certain coefficients."""
    d = round((b**2 - 4 * a * c)**(1.0 / 2), 2)
    return -b / (2 * a) if d == 0 else ((-2 * c) / (b + d), (-2 * c) / (b - d))


def relative_error(desired, obtained):
    """Presents the relative error according to results
    obtained from floating point manipulation, and a real value."""
    return abs((desired - obtained) / obtained)


origx1, origx2 = -0.0161072, -62.0839

quadx1, quadx2 = quad_eq(1, 62.1, 1)
print(round(quadx1, 4), round(quadx2, 2))
print(relative_error(origx1, round(quadx1, 4)),
      relative_error(origx2, round(quadx2, 2)))

ratx1, ratx2 = rat_quad_eq(1, 62.1, 1)
print(round(ratx1, 4), round(ratx2, 2))
print(relative_error(origx1, round(ratx1, 4)),
      relative_error(origx2, round(ratx2, 2)))
