use crate::common::point::Point;
use nalgebra::{Matrix3, Matrix3x1};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn neumann(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    u1: f64,
    u2: f64,
    u3: f64,
    edge_index: usize,
) -> Matrix3x1<f64> {
    let sqrt_2 = (2.0 as f64).sqrt();
    let foo = match edge_index {
        /* 0-edge (0,0) -> (1,0) */
        0 => Matrix3::new(
            -(p1.x / 3.0 + p2.x / 6.0 + p1.y / 3.0 + p2.y / 6.0) + 0.5,
            -(p1.x / 6.0 + p2.x / 3.0 + p1.y / 6.0 + p2.y / 3.0) + 0.5,
            0.0,
            p1.x / 3.0 + p2.x / 6.0,
            p1.x / 6.0 + p2.x / 3.0,
            0.0,
            p1.y / 3.0 + p2.y / 6.0,
            p1.y / 6.0 + p2.y / 3.0,
            0.0,
        ),
        /* 1-edge (1,0) -> (0,1) */
        1 => Matrix3::new(
            0.0,
            (-2.0 * p2.x - p3.x - 2.0 * p2.y - p3.y + 3.0) * sqrt_2 / 6.0,
            (-p2.x - 2.0 * p3.x - p2.y - 2.0 * p3.y + 3.0) * sqrt_2 / 6.0,
            0.0,
            (2.0 * p2.x + p3.x) * sqrt_2 / 6.0,
            (p2.x + 2.0 * p3.x) * sqrt_2 / 6.0,
            0.0,
            (2.0 * p2.y + p3.y) * sqrt_2 / 6.0,
            (p2.y + 2.0 * p3.y) * sqrt_2 / 6.0,
        ),
        /* 0-edge (0,1) -> (0,0) */
        2 => Matrix3::new(
            (p1.x / 3.0 + p3.x / 6.0 + p1.y / 3.0 + p3.y / 6.0) - 0.5,
            0.0,
            (p1.x / 6.0 + p3.x / 3.0 + p1.y / 6.0 + p3.y / 3.0) - 0.5,
            -(p1.x / 3.0 + p3.x / 6.0),
            0.0,
            -(p1.x / 6.0 + p3.x / 3.0),
            -(p1.y / 3.0 + p3.y / 6.0),
            0.0,
            -(p1.y / 6.0 + p3.y / 3.0),
        ),
        _ => panic!("Not expecting edge"),
    };
    return foo * Matrix3x1::new(u1, u2, u3);
}

#[cfg(test)]
mod neumann {
    use super::*;

    #[test]
    fn sample_1() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 0.0);
        let p3 = Point::new(1.0, 1.0);

        let u1 = 1.0;
        let u2 = 1.0;
        let u3 = 1.0;

        let n1 = neumann(&p1, &p2, &p3, u1, u2, u3, 0);
        let n2 = neumann(&p1, &p2, &p3, u1, u2, u3, 1);
        let n3 = neumann(&p1, &p2, &p3, u1, u2, u3, 2);
        println!("{}", n1);
        println!("{}", n2);
        println!("{}", n3);

        assert!(float_cmp::approx_eq!(
            f64,
            n1[(0, 0)],
            0.5,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n1[(1, 0)],
            0.5,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n1[(2, 0)],
            0.0,
            epsilon = 1.0E-14f64
        ));

        assert!(float_cmp::approx_eq!(
            f64,
            n2[(0, 0)],
            -0.7071067811865476,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n2[(1, 0)],
            1.4142135623730951,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n2[(2, 0)],
            0.7071067811865476,
            epsilon = 1.0E-14f64
        ));

        assert!(float_cmp::approx_eq!(
            f64,
            n3[(0, 0)],
            -0.00000000000000005551115123125783,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n3[(1, 0)],
            -0.5,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            n3[(2, 0)],
            -0.5,
            epsilon = 1.0E-14f64
        ));
    }
}
