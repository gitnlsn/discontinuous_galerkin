from sympy import *

x, y = symbols('x y')

xy = Matrix([
    [1],
    [x],
    [y],
])

a = Matrix([
    [1, 1, 1],
    [0, 1, 0],
    [0, 0, 1],
])

p0x, p1x, p2x, p0y, p1y, p2y = symbols('p0x p1x p2x p0y p1y p2y')

p = Matrix([
    [1, 1, 1],
    [p0x, p1x, p2x],
    [p0y, p1y, p2y]
])

foo = p * a.inv()

pprint(a.inv())
pprint(foo)
