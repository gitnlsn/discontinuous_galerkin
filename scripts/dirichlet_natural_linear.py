from sympy import *

x1, x2, x3 = symbols('x[0:3]')
y1, y2, y3 = symbols('y[0:3]')


def base_transformation():
    return Matrix([
        [1, -1, -1],
        [0, 1, 0],
        [0, 0, 1],
    ])


def coordinate_transformation(p1, p2, p3):
    return Matrix([
        [1, 0, 0],
        [p1[0], p2[0] - p1[0], p3[0] - p1[0]],
        [p1[1], p2[1] - p1[1], p3[1] - p1[1]],
    ])


def field_transformation(p1, p2, p3):
    return Matrix([
        [1, 1, 1],
        [p1[0], p2[0], p3[0]],
        [p1[1], p2[1], p3[1]],
    ])


def grad_n(n):
    grad = Matrix([
        [0, 0],
        [1, 0],
        [0, 1],
    ])

    if n == 0:
        return grad * Matrix([[0, -1]]).transpose()
    if n == 1:
        return grad * Matrix([[1, 1]]).transpose() / sqrt(2)
    if n == 2:
        return grad * Matrix([[-1, 0]]).transpose()


def elementaryIntegrals(n):
    x, y, l = symbols('x y lambda')

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
        return integrate(grad_n(n) * f1.transpose(), (x, 0, 1))
    if n == 1:
        return integrate(grad_n(n) * f2.transpose() * sqrt(2), (l, 0, 1))
    if n == 2:
        return integrate(grad_n(n) * f3.transpose(), (y, 1, 0))


p1 = Matrix([[x1, y1]])
p2 = Matrix([[x2, y2]])
p3 = Matrix([[x3, y3]])


def boundary_matrix(p1, p2, p3, n):
    return simplify(
        base_transformation()
        * coordinate_transformation(p1, p2, p3)
        * (
            elementaryIntegrals(n)
        )
        * coordinate_transformation(p1, p2, p3).transpose()
        * field_transformation(p1, p2, p3).transpose().inv()
    )


pprint(
    boundary_matrix(p1, p2, p3, 0)
    .subs(x1, 0).subs(y1, 0)
    .subs(x2, 1).subs(y2, 0)
    .subs(x3, 0).subs(y3, 1)
)

pprint(
    boundary_matrix(p1, p2, p3, 1)
    .subs(x1, 0).subs(y1, 0)
    .subs(x2, 1).subs(y2, 0)
    .subs(x3, 0).subs(y3, 1)
)

pprint(
    boundary_matrix(p1, p2, p3, 2)
    .subs(x1, 0).subs(y1, 0)
    .subs(x2, 1).subs(y2, 0)
    .subs(x3, 0).subs(y3, 1)
)
