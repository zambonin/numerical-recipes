from functools import reduce

for n in [2, 3, 10, 16]:
    values = reduce(lambda x, _: x + [(n + 1) * x[-1] - 1], range(10), [1 / n])
    print("\n".join("{:.30f}".format(i) for i in values))
