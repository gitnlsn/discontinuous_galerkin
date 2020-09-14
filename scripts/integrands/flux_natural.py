from sympy import *

x, y, l = symbols('x y lambda')


def elementaryIntegrals(n):
    f1 = Matrix([
        [1],
        [x],
        [0],
    ])

    f2 = Matrix([
        [1],
        [1-l],
        [l],
    ])

    f3 = Matrix([
        [1],
        [0],
        [y],
    ])

    if n == 0:
        return integrate(f1, (x, 0, 1))
    if n == 1:
        return integrate(f2 * sqrt(2), (l, 0, 1))
    if n == 2:
        return integrate(f3, (y, 1, 0))


def coordinate_transformation(p1, p2, p3):
    return Matrix([
        [1, 0, 0],
        [p1[0], p2[0] - p1[0], p3[0] - p1[0]],
        [p1[1], p2[1] - p1[1], p3[1] - p1[1]],
    ])


def base_transformation():
    return Matrix([
        [1, -1, -1],
        [0, 1, 0],
        [0, 0, 1],
    ])


def field_transformation(p1, p2, p3):
    return Matrix([
        [1, 1, 1],
        [p1[0], p2[0], p3[0]],
        [p1[1], p2[1], p3[1]],
    ])


def grad(n):
    if n == 0:
        return Matrix([[0, -1]]) * Matrix([
            [0, 0],
            [1, 0],
            [0, 1],
        ]).transpose()
    if n == 1:
        return Matrix([[1, 1]]) / sqrt(2) * Matrix([
            [0, 0],
            [1, 0],
            [0, 1],
        ]).transpose()
    if n == 2:
        return Matrix([[-1, 0]]) * Matrix([
            [0, 0],
            [1, 0],
            [0, 1],
        ]).transpose()


def flux_matrix(p1, p2, p3, n):
    return simplify(
        base_transformation()
        * coordinate_transformation(p1, p2, p3)
        * elementaryIntegrals(n)
        * grad(n)
    )


x1, x2, x3, x4, x5, x6 = symbols('x[0:6]')
y1, y2, y3, y4, y5, y6 = symbols('y[0:6]')

p1 = Matrix([[x1, y1]])
p2 = Matrix([[x2, y2]])
p3 = Matrix([[x3, y3]])
p4 = Matrix([[x4, y4]])
p5 = Matrix([[x5, y5]])
p6 = Matrix([[x6, y6]])

pprint(flux_matrix(p1, p2, p3, 0))
pprint(flux_matrix(p1, p2, p3, 1))
pprint(flux_matrix(p1, p2, p3, 2))
