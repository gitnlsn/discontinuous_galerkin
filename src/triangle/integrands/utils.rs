use crate::common::point::Point;
use nalgebra::{Matrix1x2, Matrix3, Matrix3x1, Matrix3x2};

/**
 * Calculates integrals:
 *  [1, x, 0] {x, 0, 1}
 *  [1, 1-l, l] * sqrt(2) {l, 0, 1}
 *  [1, 0, y] {y, 1, 0}
 */
pub fn elementary_intergral(n: usize) -> Matrix3x1<f64> {
    match n {
        0 => Matrix3x1::new(1.0, 0.5, 0.0),
        1 => Matrix3x1::new(
            (2.0 as f64).sqrt(),
            (2.0 as f64).sqrt() / 2.0,
            (2.0 as f64).sqrt() / 2.0,
        ),
        2 => Matrix3x1::new(1.0, 0.0, -0.5),
        _ => panic!("Not expected integral definition greater than 2"),
    }
}

/**
 * Maps [1,x,y] into [1-x-y, x, y]
 */
pub fn base_transformation() -> Matrix3<f64> {
    Matrix3::new(1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
}

/**
 * Projects u(x,y) as sum of [1,x,y] over [u0, u1, u2]
 *  u(x,y) = [1] * [1, x0, y0] ^ -1 * [u0]
 *           [x]   [1, x1, y1]        [u1]
 *           [y]   [1, x2, y2]        [u2]
 */
pub fn field_transformation(p1: &Point, p2: &Point, p3: &Point) -> Matrix3<f64> {
    Matrix3::new(
        1.0, p1.x, p1.y,
        1.0, p2.x, p2.y,
        1.0, p3.x, p3.y,
    )
}

/**
 * Coordinate transformation for default triangle ((0,0),(1,0), (0,1))
 */
pub fn coordinante_transformation(p1: &Point, p2: &Point, p3: &Point) -> Matrix3<f64> {
    Matrix3::new(
        1.0,
        0.0,
        0.0,
        p1.x,
        p2.x - p1.x,
        p3.x - p1.x,
        p1.y,
        p2.y - p1.y,
        p3.y - p1.y,
    )
}

/**
 * Gradient has constant direction for default triangle ((0,0),(1,0), (0,1))
 */
pub fn grad() -> Matrix3x2<f64> {
    Matrix3x2::new(-1.0, -1.0, 1.0, 0.0, 0.0, 1.0)
}

/**
 * Outer normal versor at boundaries for default triangle ((0,0),(1,0), (0,1))
 */
pub fn normal(n: usize) -> Matrix1x2<f64> {
    match n {
        0 => Matrix1x2::new(0.0, -1.0),
        1 => Matrix1x2::new(1.0 / (2.0 as f64).sqrt(), 1.0 / (2.0 as f64).sqrt()),
        2 => Matrix1x2::new(-1.0, 0.0),
        _ => panic!("Not expected normal with index greater than 2"),
    }
}