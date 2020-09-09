use crate::common::point::Point;
use nalgebra::Matrix3;

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn penalty(p1: &Point, p2: &Point) -> Matrix3<f64> {
    Matrix3::new(
        -p1.x / 3.0 - p2.x / 6.0 - p1.y / 3.0 - p2.y / 6.0 + 0.5,
        -p1.x / 6.0 - p2.x / 3.0 - p1.y / 6.0 - p2.y / 3.0 + 0.5,
        0.0,
        p1.x / 3.0 + p2.x / 6.0,
        p1.x / 6.0 + p2.x / 3.0,
        0.0,
        p1.y / 3.0 + p2.y / 6.0,
        p1.y / 6.0 + p2.y / 3.0,
        0.0,
    )
}

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn half_penalty(p1: &Point, p2: &Point) -> Matrix3<f64> {
    Matrix3::new(
        -p1.x / 6.0 - p2.x / 12.0 - p1.y / 6.0 - p2.y / 12.0 + 0.25,
        -p1.x / 12.0 - p2.x / 6.0 - p1.y / 12.0 - p2.y / 6.0 + 0.25,
        0.0,
        p1.x / 6.0 + p2.x / 12.0,
        p1.x / 12.0 + p2.x / 6.0,
        0.0,
        p1.y / 6.0 + p2.y / 12.0,
        p1.y / 12.0 + p2.y / 6.0,
        0.0,
    )
}
