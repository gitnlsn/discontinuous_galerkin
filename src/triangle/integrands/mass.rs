extern crate nalgebra;

use crate::common::point::Point;

use nalgebra::{Matrix2, Matrix2x1, Matrix3};

pub fn matrix(p1: &Point, p2: &Point, p3: &Point) -> Matrix3<f64> {
    let jacobian = jacobian(&p1, &p2, &p3);
    return Matrix3::new(
        mass_ij(&jacobian, 0, 0),
        mass_ij(&jacobian, 0, 1),
        mass_ij(&jacobian, 0, 2),
        mass_ij(&jacobian, 1, 0),
        mass_ij(&jacobian, 1, 1),
        mass_ij(&jacobian, 1, 2),
        mass_ij(&jacobian, 2, 0),
        mass_ij(&jacobian, 2, 1),
        mass_ij(&jacobian, 2, 2),
    );
}

fn jacobian(p1: &Point, p2: &Point, p3: &Point) -> Matrix2<f64> {
    return Matrix2::new(p2.x - p1.x, p3.x - p1.x, p2.y - p1.y, p3.y - p1.y);
}

fn dphi(k: usize) -> Matrix2x1<f64> {
    /*
        Warning(!): hardcoded derivatives
    */
    match k {
        0 /* phi(x,y) = 1-x-y */ => Matrix2x1::new( -1.0, -1.0),
        1 /* phi(x,y) = x     */ => Matrix2x1::new( 1.0, 0.0),
        2 /* phi(x,y) = y     */ => Matrix2x1::new( 0.0, 1.0),
        _ => panic!("Not expected to request later polynomial"),
    }
}

fn mass_ij(jaco: &Matrix2<f64>, i: usize, j: usize) -> f64 {
    let jaco_inv = jaco.try_inverse().unwrap();

    let dphi_i = dphi(i);
    let dphi_j = dphi(j);

    let di_dx = dphi_i[0];
    let di_dy = dphi_i[1];
    let dj_dx = dphi_j[0];
    let dj_dy = dphi_j[1];

    let j11 = jaco_inv[(0, 0)]; /* 1,1 */
    let j12 = jaco_inv[(0, 1)]; /* 1,2 */
    let j21 = jaco_inv[(1, 0)]; /* 2,1 */
    let j22 = jaco_inv[(1, 1)]; /* 2,2 */

    let coef = j11 * j11 * di_dx * dj_dx
        + j11 * j12 * (di_dy * dj_dx + di_dx * dj_dy)
        + j12 * j12 * di_dy * dj_dy
        + j21 * j21 * di_dx * dj_dx
        + j21 * j22 * (di_dy * dj_dx + di_dx * dj_dy)
        + j22 * j22 * di_dy * dj_dy;

    return coef * jaco.determinant();
}

#[cfg(test)]
mod mass_matrix {
    use super::*;

    #[test]
    fn sample_1() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let p3 = Point::new(0.0, 1.0);

        let jaco = jacobian(&p1, &p2, &p3);

        assert_eq!(mass_ij(&jaco, 0, 0), 2.0);
        assert_eq!(mass_ij(&jaco, 0, 1), -1.0);
        assert_eq!(mass_ij(&jaco, 0, 2), -1.0);
        assert_eq!(mass_ij(&jaco, 1, 0), -1.0);
        assert_eq!(mass_ij(&jaco, 1, 1), 1.0);
        assert_eq!(mass_ij(&jaco, 1, 2), 0.0);
        assert_eq!(mass_ij(&jaco, 2, 0), -1.0);
        assert_eq!(mass_ij(&jaco, 2, 1), 0.0);
        assert_eq!(mass_ij(&jaco, 2, 2), 1.0);
    }

    #[test]
    fn sample_2() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(0.1, 0.0);
        let p3 = Point::new(0.0, 0.1);

        let jaco = jacobian(&p1, &p2, &p3);

        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 0, 0),
            2.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 0, 1),
            -1.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 0, 2),
            -1.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 1, 0),
            -1.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 1, 1),
            1.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 1, 2),
            0.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 2, 0),
            -1.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 2, 1),
            0.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            mass_ij(&jaco, 2, 2),
            1.0,
            epsilon = 1.0E-14f64
        ));
    }
}
