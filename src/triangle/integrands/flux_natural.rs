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
            0.0,
            0.0,
            (p1.x + p2.x + p1.y + p2.y) / 2.0 - 1.0,
            0.0,
            0.0,
            -(p1.x + p2.x) / 2.0,
            0.0,
            0.0,
            -(p1.y + p2.y) / 2.0,
        ),
        1 => Matrix3::new(
            0.0,
            -(p2.x + p3.x + p2.y + p3.y) / 2.0 + 1.0,
            -(p2.x + p3.x + p2.y + p3.y) / 2.0 + 1.0,
            0.0,
            (p2.x + p3.x) / 2.0,
            (p2.x + p3.x) / 2.0,
            0.0,
            (p2.y + p3.y) / 2.0,
            (p2.y + p3.y) / 2.0,
        ),
        2 => Matrix3::new(
            0.0,
            -(p1.x + p3.x + p1.y + p3.y) / 2.0 + 1.0,
            0.0,
            0.0,
            (p1.x + p3.x) / 2.0,
            0.0,
            0.0,
            (p1.y + p3.y) / 2.0,
            0.0,
        ),
        _ => panic!("Not expected to have edge index greater than 2 for triangle"),
    };
    return foo / 2.0
        * utils::coordinante_transformation(p1, p2, p3).transpose()
        * utils::field_transformation(p4, p5, p6)
            .try_inverse()
            .unwrap();
}

#[cfg(test)]
mod flux_matrix {
    use super::*;

    #[test]
    fn sample_1() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let p3 = Point::new(0.0, 1.0);

        let edge_1 = half_flux(&p1, &p2, &p3, &p1, &p2, &p3, 0);
        assert_eq!(edge_1[(0, 0)], 0.25);
        assert_eq!(edge_1[(0, 1)], 0.0);
        assert_eq!(edge_1[(0, 2)], -0.25);

        assert_eq!(edge_1[(1, 0)], 0.25);
        assert_eq!(edge_1[(1, 1)], 0.0);
        assert_eq!(edge_1[(1, 2)], -0.25);

        assert_eq!(edge_1[(2, 0)], 0.0);
        assert_eq!(edge_1[(2, 1)], 0.0);
        assert_eq!(edge_1[(2, 2)], 0.0);

        let edge_2 = half_flux(&p1, &p2, &p3, &p1, &p2, &p3, 1);
        assert_eq!(edge_2[(0, 0)], 0.0);
        assert_eq!(edge_2[(0, 1)], 0.0);
        assert_eq!(edge_2[(0, 2)], 0.0);

        assert_eq!(edge_2[(1, 0)], -0.5);
        assert_eq!(edge_2[(1, 1)], 0.25);
        assert_eq!(edge_2[(1, 2)], 0.25);

        assert_eq!(edge_2[(2, 0)], -0.5);
        assert_eq!(edge_2[(2, 1)], 0.25);
        assert_eq!(edge_2[(2, 2)], 0.25);

        let edge_3 = half_flux(&p1, &p2, &p3, &p1, &p2, &p3, 2);
        assert_eq!(edge_3[(0, 0)], -0.25);
        assert_eq!(edge_3[(0, 1)], 0.25);
        assert_eq!(edge_3[(0, 2)], 0.0);

        assert_eq!(edge_3[(1, 0)], 0.0);
        assert_eq!(edge_3[(1, 1)], 0.0);
        assert_eq!(edge_3[(1, 2)], 0.0);

        assert_eq!(edge_3[(2, 0)], -0.25);
        assert_eq!(edge_3[(2, 1)], 0.25);
        assert_eq!(edge_3[(2, 2)], 0.0);
    }
}
