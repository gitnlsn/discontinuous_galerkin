use crate::common::point::Point;
use nalgebra::{Matrix3};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn flux(p1: &Point, p2: &Point, p3: &Point) -> Matrix3<f64> {
    Matrix3::new(
        (-p1.x + p3.x - p1.y + p3.y) / 2.0,
        (-p1.x + p3.x - p1.y + p3.y) / 2.0,
        0.0,
        (p1.x - p3.x) / 2.0,
        (p1.x - p3.x) / 2.0,
        0.0,
        (p1.y - p3.y) / 2.0,
        (p1.y - p3.y) / 2.0,
        0.0,
    )
}


/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn half_flux(p1: &Point, p2: &Point, p3: &Point) -> Matrix3<f64> {
    Matrix3::new(
        (-p1.x + p3.x - p1.y + p3.y) / 4.0,
        (-p1.x + p3.x - p1.y + p3.y) / 4.0,
        0.0,
        (p1.x - p3.x) / 4.0,
        (p1.x - p3.x) / 4.0,
        0.0,
        (p1.y - p3.y) / 4.0,
        (p1.y - p3.y) / 4.0,
        0.0,
    )
}
