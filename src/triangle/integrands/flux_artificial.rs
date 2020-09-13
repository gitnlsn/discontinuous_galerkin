use crate::common::point::Point;
use crate::triangle::integrands::utils;
use nalgebra::Matrix3;

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn half_flux(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    p4: &Point,
    p5: &Point,
    p6: &Point,
    edge_index: usize,
) -> Matrix3<f64> {
    let foo = match edge_index {
        0 => Matrix3::new(
            -p1.x + p3.x - p1.y + p3.y,
            (-p1.x + p3.x - p1.y + p3.y) / 2.0,
            0.0,
            p1.x - p3.x,
            (p1.x - p3.x) / 2.0,
            0.0,
            p1.y - p3.y,
            (p1.y - p3.y) / 2.0,
            0.0,
        ),
        1 => Matrix3::new(
            2.0 * p1.x - p2.x - p3.x + 2.0 * p1.y - p2.y - p3.y,
            (2.0 * p1.x - p2.x - p3.x + 2.0 * p1.y - p2.y - p3.y) / 2.0,
            (2.0 * p1.x - p2.x - p3.x + 2.0 * p1.y - p2.y - p3.y) / 2.0,
            -2.0 * p1.x + p2.x + p3.x,
            (-2.0 * p1.x + p2.x + p3.x) / 2.0,
            (-2.0 * p1.x + p2.x + p3.x) / 2.0,
            -2.0 * p1.y + p2.y + p3.y,
            (-2.0 * p1.y + p2.y + p3.y) / 2.0,
            (-2.0 * p1.y + p2.y + p3.y) / 2.0,
        ),
        2 => Matrix3::new(
            p1.x - p2.x + p1.y - p2.y,
            0.0,
            (p1.x - p2.x + p1.y - p2.y) / 2.0,
            -p1.x + p2.x,
            0.0,
            (-p1.x + p2.x) / 2.0,
            -p1.y + p2.y,
            0.0,
            (-p1.y + p2.y) / 2.0,
        ),
        _ => panic!("Not expected edge index greater than 2 for triangle"),
    };
    return foo / 2.0
        * utils::coordinante_transformation(p1, p2, p3).transpose()
        * utils::field_transformation(p4, p5, p6)
            .try_inverse()
            .unwrap();
}
