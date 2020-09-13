use crate::common::point::Point;
use nalgebra::{Matrix3, Matrix3x1};

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn dirichlet_linear_natural(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    u1: f64,
    u2: f64,
    u3: f64,
    edge_index: usize,
) -> Matrix3x1<f64> {
    let a = match edge_index {
        /* 0-edge (0,0) -> (1,0) */
        0 => Matrix3::new(
            (-p1.x + p3.x - p1.y + p3.y) / 2.0,
            (-p1.x + p3.x - p1.y + p3.y) / 2.0,
            0.0,
            (p1.x - p3.x) / 2.0,
            (p1.x - p3.x) / 2.0,
            0.0,
            (p1.y - p3.y) / 2.0,
            (p1.y - p3.y) / 2.0,
            0.0,
        ),
        /* 1-edge (1,0) -> (0,1) */
        1 => Matrix3::new(
            0.0,
            p1.x - p2.x / 2.0 - p3.x / 2.0 + p1.y - p2.y / 2.0 - p3.y / 2.0,
            p1.x - p2.x / 2.0 - p3.x / 2.0 + p1.y - p2.y / 2.0 - p3.y / 2.0,
            0.0,
            -p1.x + p2.x / 2.0 + p3.x / 2.0,
            -p1.x + p2.x / 2.0 + p3.x / 2.0,
            0.0,
            -p1.y + p2.y / 2.0 + p3.y / 2.0,
            -p1.y + p2.y / 2.0 + p3.y / 2.0,
        ),
        /* 0-edge (0,1) -> (0,0) */
        2 => Matrix3::new(
            (p1.x - p2.x + p1.y - p2.y) / 2.0,
            0.0,
            (p1.x - p2.x + p1.y - p2.y) / 2.0,
            (-p1.x + p2.x) / 2.0,
            0.0,
            (-p1.x + p2.x) / 2.0,
            (-p1.y + p2.y) / 2.0,
            0.0,
            (-p1.y + p2.y) / 2.0,
        ),
        _ => panic!("Not expecting edge"),
    };
    return a * Matrix3x1::new(u1, u2, u3);
}

pub fn dirichlet_bilinear_penalty(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    edge_index: usize,
) -> Matrix3<f64> {
    let sqrt_2 = (2.0 as f64).sqrt();
    match edge_index {
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
    }
}

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn dirichlet_linear_penalty(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    u1: f64,
    u2: f64,
    u3: f64,
    edge_index: usize,
) -> Matrix3x1<f64> {
    dirichlet_bilinear_penalty(p1, p2, p3, edge_index) * Matrix3x1::new(u1, u2, u3)
}

#[cfg(test)]
mod linear_natural {
    use super::*;
    use nalgebra::DMatrix;
    use std::rc::Rc;

    #[test]
    fn sample_1() {
        /*
            Sigma-0 Penalty matrices will lead to strongly 
            imposed dirichlet conditions in single element
        */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(0.0, 1.0));

        /* Assembling system */
        let mut system_matrix = DMatrix::<f64>::zeros(3, 3);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 0);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 1);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 2);

        /* Extern matrix */
        let mut extern_matrix = DMatrix::<f64>::zeros(3, 1);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 0);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 1);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 2);

        let solution = system_matrix.try_inverse().unwrap() * extern_matrix;
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(0, 0)],
            2.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(1, 0)],
            0.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(2, 0)],
            0.0,
            epsilon = 1.0E-14f64
        ));
    }

    #[test]
    fn sample_2() {
        /*
            Sigma-0 Penalty matrices will lead to strongly 
            imposed dirichlet conditions in single element

            Independently of the coodinates of the element
        */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(3.3, 0.0));
        let p3 = Rc::new(Point::new(0.0, 8.43214));

        /* Assembling system */
        let mut system_matrix = DMatrix::<f64>::zeros(3, 3);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 0);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 1);
        system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 2);

        /* Extern matrix */
        let mut extern_matrix = DMatrix::<f64>::zeros(3, 1);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 0);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 1);
        extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 2);

        let solution = system_matrix.try_inverse().unwrap() * extern_matrix;
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(0, 0)],
            2.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(1, 0)],
            0.0,
            epsilon = 1.0E-14f64
        ));
        assert!(float_cmp::approx_eq!(
            f64,
            solution[(2, 0)],
            0.0,
            epsilon = 1.0E-14f64
        ));
    }

    #[test]
    fn sample_3() {
        /*
            Solution will converge if sigma increases
        */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(0.0, 1.0));

        fn solve(p1: &Rc<Point>, p2: &Rc<Point>, p3: &Rc<Point>, sigma: f64) -> DMatrix<f64> {
            /* Assembling system */
            let mut system_matrix = DMatrix::<f64>::zeros(3, 3);
            system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 0) * sigma;
            system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 1) * sigma;
            system_matrix += dirichlet_bilinear_penalty(&p1, &p2, &p3, 2) * sigma;

            /* Extern matrix */
            let mut extern_matrix = DMatrix::<f64>::zeros(3, 1);
            extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 0) * sigma;
            extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 1) * sigma;
            extern_matrix += dirichlet_linear_penalty(&p1, &p2, &p3, 2.0, 0.0, 0.0, 2) * sigma;

            extern_matrix += dirichlet_linear_natural(&p1, &p2, &p3, 2.0, 0.0, 0.0, 0);
            extern_matrix += dirichlet_linear_natural(&p1, &p2, &p3, 2.0, 0.0, 0.0, 1);
            extern_matrix += dirichlet_linear_natural(&p1, &p2, &p3, 2.0, 0.0, 0.0, 2);

            return system_matrix.try_inverse().unwrap() * extern_matrix;
        }

        let sol_1 = solve(&p1, &p2, &p3, 100.0);
        let sol_2 = solve(&p1, &p2, &p3, 1000.0);
        let sol_3 = solve(&p1, &p2, &p3, 10000.0);
        let sol_4 = solve(&p1, &p2, &p3, 100000.0);

        let delta_1 = ((&sol_2 - &sol_1).transpose() * (&sol_2 - &sol_1))[(0, 0)];
        let delta_2 = ((&sol_3 - &sol_2).transpose() * (&sol_3 - &sol_2))[(0, 0)];
        let delta_3 = ((&sol_4 - &sol_3).transpose() * (&sol_4 - &sol_3))[(0, 0)];
        
        assert!(delta_2 < delta_1);
        assert!(delta_3 < delta_2);
    }
}
