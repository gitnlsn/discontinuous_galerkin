from sympy import *

x, y = symbols('x y')


def phi(n):
    if n == 0:
        return 1 - x - y
    elif n == 1:
        return x
    elif n == 2:
        return y
    else:
        pass


def jacobian(p1, p2, p3):
    return Matrix([
        [p2[0] - p1[0], p3[0] - p1[0]],
        [p2[1] - p1[1], p3[1] - p1[1]],
    ])


def mass_ij(i, j, p1, p2, p3):
    phi_i = phi(i)
    phi_j = phi(j)

    di_dx = diff(phi_i, x)
    di_dy = diff(phi_i, y)
    dj_dx = diff(phi_j, x)
    dj_dy = diff(phi_j, y)

    J = jacobian(p1, p2, p3)
    J_inv = J.inv()
    j11 = J_inv.row(0).col(0)
    j12 = J_inv.row(0).col(1)
    j21 = J_inv.row(1).col(0)
    j22 = J_inv.row(1).col(1)
    j_det = J.det()

    coef = j11 * j11 * di_dx * dj_dx                    \
        + j11 * j12 * (di_dy * dj_dx + di_dx * dj_dy)   \
        + j12 * j12 * di_dy * dj_dy                     \
        + j21 * j21 * di_dx * dj_dx                     \
        + j21 * j22 * (di_dy * dj_dx + di_dx * dj_dy)   \
        + j22 * j22 * di_dy * dj_dy

    return integrate(integrate(coef * j_det, (x, 0, 1)), (y, 0, 1))


def mass(p1, p2, p3):
    return Matrix([
        [mass_ij(0, 0, p1, p2, p3), mass_ij(
            0, 1, p1, p2, p3), mass_ij(0, 2, p1, p2, p3)],
        [mass_ij(1, 0, p1, p2, p3), mass_ij(
            1, 1, p1, p2, p3), mass_ij(1, 2, p1, p2, p3)],
        [mass_ij(2, 0, p1, p2, p3), mass_ij(
            2, 1, p1, p2, p3), mass_ij(2, 2, p1, p2, p3)],
    ])


p0 = Matrix([[0, 0]])
p1 = Matrix([[1, 0]])
p2 = Matrix([[0, 1]])

pprint(mass(p0, p1, p2))
