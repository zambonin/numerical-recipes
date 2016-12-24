from numpy import float32 as single

h = 1 / 2
x = 2 / 3 - h
y = 3 / 5 - h
e = 3 * x - h
f = 5 * y - h
g = e / f

for i in [h, x, y, e, f, g]:
    print("{:.55f}\n{:.55f}\n".format(single(i), i))
