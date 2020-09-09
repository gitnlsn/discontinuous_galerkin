use crate::common::point::Point;
use nalgebra::{Matrix3x1};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn neumann(p1: &Point, p2: &Point, gn: f64) -> Matrix3x1<f64> {
    Matrix3x1::new(
        -p1.x / 2.0 - p2.x / 2.0 - p1.y / 2.0 - p2.y / 2.0 + 0.5,
        p1.x / 2.0 + p2.x / 2.0,
        p1.y / 2.0 + p2.y / 2.0,
    ) * gn
}
