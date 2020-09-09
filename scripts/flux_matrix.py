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


def grad(p1, p2, p3):
    return Matrix([
        [0, 0],
        [1, 0],
        [0, 1],
    ]).transpose()                                          \
        * coordinate_transformation(p1, p2, p3).transpose() \
        * field_transformation(p1, p2, p3).transpose().inv()


x1, y1, x2, y2, x3, y3, x4, y4 = symbols('x1,y1, x2,y2, x3,y3, x4,y4')

p1 = Matrix([[x1, y1]])
p2 = Matrix([[x2, y2]])
p3 = Matrix([[x3, y3]])
p4 = Matrix([[x4, y4]])

f_matrix = simplify(
    base_transformation()
    * coordinate_transformation(p2, p4, p1)
    * elementaryIntegrals(0)
    * Matrix([[0, -1]])
    * grad(p2, p4, p1)
)

f_matrix_opp = simplify(
    base_transformation()
    * coordinate_transformation(p4, p2, p3)
    * elementaryIntegrals(0)
    * Matrix([[0, -1]])
    * grad(p2, p4, p1)
)

pprint(f_matrix)
pprint(f_matrix_opp)
pprint(
    f_matrix
    .subs(x1, 0).subs(y1, 0)
    .subs(x2, 1).subs(y2, 0)
    .subs(x3, 1).subs(y3, 1)
    .subs(x4, 0).subs(y4, 1)
)
