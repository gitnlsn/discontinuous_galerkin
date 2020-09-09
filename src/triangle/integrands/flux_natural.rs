use crate::common::point::Point;
use nalgebra::{Matrix1x2, Matrix3, Matrix3x1, Matrix3x2};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn flux(p1: &Point, p2: &Point) -> Matrix3<f64> {
    Matrix3::new(
        -(p1.x + p2.x + p1.y + p2.y) / 2.0 + 1.0,
        0.0,
        (p1.x + p2.x + p1.y + p2.y) / 2.0 - 1.0,
        (p1.x + p2.x) / 2.0,
        0.0,
        -(p1.x + p2.x) / 2.0,
        (p1.y + p2.y) / 2.0,
        0.0,
        -(p1.y + p2.y) / 2.0,
    )
}

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn half_flux(p1: &Point, p2: &Point) -> Matrix3<f64> {
    Matrix3::new(
        -(p1.x + p2.x + p1.y + p2.y) / 4.0 + 0.5,
        0.0,
        (p1.x + p2.x + p1.y + p2.y) / 4.0 - 0.5,
        (p1.x + p2.x) / 4.0,
        0.0,
        -(p1.x + p2.x) / 4.0,
        (p1.y + p2.y) / 4.0,
        0.0,
        -(p1.y + p2.y) / 4.0,
    )
}

#[cfg(test)]
mod flux_matrix {
    use super::*;

    #[test]
    fn sample_1() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let f_matrix = flux(&p1, &p2);

        assert_eq!(f_matrix[(0, 0)], 0.5);
        assert_eq!(f_matrix[(0, 1)], 0.0);
        assert_eq!(f_matrix[(0, 2)], -0.5);

        assert_eq!(f_matrix[(1, 0)], 0.5);
        assert_eq!(f_matrix[(1, 1)], 0.0);
        assert_eq!(f_matrix[(1, 2)], -0.5);

        assert_eq!(f_matrix[(2, 0)], 0.0);
        assert_eq!(f_matrix[(2, 1)], 0.0);
        assert_eq!(f_matrix[(2, 2)], 0.0);
    }

    #[test]
    fn sample_2() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let f_matrix = half_flux(&p1, &p2);

        assert_eq!(f_matrix[(0, 0)], 0.25);
        assert_eq!(f_matrix[(0, 1)], 0.0);
        assert_eq!(f_matrix[(0, 2)], -0.25);

        assert_eq!(f_matrix[(1, 0)], 0.25);
        assert_eq!(f_matrix[(1, 1)], 0.0);
        assert_eq!(f_matrix[(1, 2)], -0.25);

        assert_eq!(f_matrix[(2, 0)], 0.0);
        assert_eq!(f_matrix[(2, 1)], 0.0);
        assert_eq!(f_matrix[(2, 2)], 0.0);
    }
}
