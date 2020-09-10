use crate::common::point::Point;
use nalgebra::{Matrix3, Matrix3x1};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn neumann(p1: &Point, p2: &Point, p3: &Point, u1: f64, u2: f64) -> Matrix3x1<f64> {
    Matrix3::new(
        -(p1.x / 3.0 + p2.x / 6.0 + p1.y / 3.0 + p2.y / 6.0) + 0.5,
        -(p1.x / 6.0 + p2.x / 3.0 + p1.y / 6.0 + p2.y / 3.0) + 0.5,
        0.0,
        p1.x / 3.0 + p2.x / 6.0,
        p1.x / 6.0 + p2.x / 3.0,
        0.0,
        p1.y / 3.0 + p2.y / 6.0,
        p1.y / 6.0 + p2.y / 3.0,
        0.0,
    ) * Matrix3x1::new(u1, u2, 0.0)
}
