use nalgebra::DMatrix;

use crate::triangle::system_builder::{
    domain::Domain,
    external::{dirichlet, neumann},
    internal::{flux_artificial, flux_natural, jump_penalty, mass},
};

pub fn build(sigma: f64, domain: &Domain) -> (DMatrix<f64>, DMatrix<f64>) {
    let system_size = domain.elements.len() * 3;

    let mut system_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(system_size, system_size);
    let mut extern_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(system_size, 1);

    mass::build(&mut system_matrix, domain).expect("Failed to build mass matrix");
    flux_natural::build(&mut system_matrix, domain).expect("Failed to set natural flux");
    flux_artificial::build(&mut system_matrix, domain)
        .expect("Failed to set artificial anti-symmetric flux");

    jump_penalty::build(&mut system_matrix, sigma, domain).expect("Failed to set inter element penalties");
    dirichlet::build(&mut system_matrix, &mut extern_matrix, sigma, domain).expect("Failed to set dirichlet constraints");
    neumann::build(&mut extern_matrix, domain).expect("Failed to set neumann constraints");

    return (system_matrix, extern_matrix);
}

#[cfg(test)]
mod system_build {
    use super::*;
    use crate::common::{edge::Edge, point::Point};
    use crate::triangle::element::TriangleElementL1;

    use std::rc::Rc;

    #[test]
    fn sample_1() {
        /* 6 triangle hexagon */
        let p1 = Rc::new(Point::new(2.0, 1.0));
        let p2 = Rc::new(Point::new(3.0, 2.0));
        let p3 = Rc::new(Point::new(3.0, 4.0));
        let p4 = Rc::new(Point::new(2.0, 5.0));
        let p5 = Rc::new(Point::new(1.0, 4.0));
        let p6 = Rc::new(Point::new(1.0, 2.0));
        let p7 = Rc::new(Point::new(2.0, 3.0));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p7));
        let t2 = Rc::new(TriangleElementL1::new(&p2, &p3, &p7));
        let t3 = Rc::new(TriangleElementL1::new(&p3, &p4, &p7));
        let t4 = Rc::new(TriangleElementL1::new(&p4, &p5, &p7));
        let t5 = Rc::new(TriangleElementL1::new(&p5, &p6, &p7));
        let t6 = Rc::new(TriangleElementL1::new(&p6, &p1, &p7));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);
        domain.insert_element(&t3);
        domain.insert_element(&t4);
        domain.insert_element(&t5);
        domain.insert_element(&t6);

        let e1 = Rc::new(Edge::new(&p1, &p2));
        let e2 = Rc::new(Edge::new(&p2, &p3));
        let e3 = Rc::new(Edge::new(&p3, &p4));
        let e4 = Rc::new(Edge::new(&p4, &p5));
        let e5 = Rc::new(Edge::new(&p5, &p6));
        let e6 = Rc::new(Edge::new(&p6, &p1));

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 1.0]);
        domain.insert_dirichlet_constraint(&e2, vec![1.0, 1.0]);
        domain.insert_dirichlet_constraint(&e3, vec![1.0, 0.0]);
        domain.insert_dirichlet_constraint(&e4, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e5, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e6, vec![0.0, 0.0]);

        fn solve(sigma: f64, domain: &Domain) -> f64 {
            let (system_matrix, extern_matrix) = build(sigma, &domain);
            let inverse = system_matrix.try_inverse().unwrap();
            let answer = inverse * extern_matrix;

            let expected = DMatrix::<f64>::from_row_slice(
                18,
                1,
                &vec![
                    0.0, 1.0, -1.0, /*  */
                    1.0, 1.0, -1.0, /*  */
                    1.0, 0.0, -1.0, /*  */
                    0.0, 0.0, -1.0, /*  */
                    0.0, 0.0, -1.0, /*  */
                    0.0, 0.0, -1.0, /*  */
                ],
            );
            let error = ((&answer - &expected).transpose() * (&answer - &expected))[(0, 0)].sqrt();
            // println!("{} {}", answer, error);
            return error;
        }
        let err1 = solve(10.0, &domain);
        let err2 = solve(100.0, &domain);
        let err3 = solve(1000.0, &domain);
        let err4 = solve(10000.0, &domain);
        let err5 = solve(100000.0, &domain);

        assert!(err2 < err1);
        assert!(err3 < err2);
        assert!(err4 < err3);
        assert!(err5 < err4);
    } /* end - hexagon */

    #[test]
    fn sample_2() {
        /* square into triangles */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let e1 = Rc::new(Edge::new(&p1, &p2));
        let e2 = Rc::new(Edge::new(&p2, &p3));
        let e3 = Rc::new(Edge::new(&p3, &p4));
        let e4 = Rc::new(Edge::new(&p4, &p1));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p4, &p2, &p3));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e2, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e3, vec![0.0, 1.0]);
        domain.insert_dirichlet_constraint(&e4, vec![1.0, 0.0]);

        fn solve(sigma: f64, domain: &Domain) -> f64 {
            let (system_matrix, extern_matrix) = build(sigma, &domain);
            let inverse = system_matrix.try_inverse().unwrap();
            let answer = inverse * extern_matrix;

            let expected =
                DMatrix::<f64>::from_row_slice(6, 1, &vec![0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
            let error = ((&answer - &expected).transpose() * (&answer - &expected))[(0, 0)].sqrt();
            // println!("{} {}", answer, error);
            return error;
        }
        let err1 = solve(10.0, &domain);
        let err2 = solve(100.0, &domain);
        let err3 = solve(1000.0, &domain);
        let err4 = solve(10000.0, &domain);
        let err5 = solve(100000.0, &domain);

        assert!(err2 < err1);
        assert!(err3 < err2);
        assert!(err4 < err3);
        assert!(err5 < err4);
    }

    #[test]
    fn sample_3() {
        /* Single triangle */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let e1 = Rc::new(Edge::new(&p1, &p2));
        let ed = Rc::new(Edge::new(&p2, &p4));
        let e4 = Rc::new(Edge::new(&p4, &p1));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 1.0]);
        domain.insert_dirichlet_constraint(&ed, vec![1.0, 0.0]);
        domain.insert_dirichlet_constraint(&e4, vec![0.0, 0.0]);

        fn solve(sigma: f64, domain: &Domain) -> f64 {
            let (system_matrix, extern_matrix) = build(sigma, &domain);
            let inverse = system_matrix.try_inverse().unwrap();
            let answer = inverse * extern_matrix;

            let expected = DMatrix::<f64>::from_row_slice(3, 1, &vec![0.0, 1.0, 0.0]);
            let error = ((&answer - &expected).transpose() * (&answer - &expected))[(0, 0)].sqrt();
            // println!("{} {}", answer, error);
            return error;
        }
        let err1 = solve(1.0, &domain);
        let err2 = solve(10.0, &domain);
        let err3 = solve(100.0, &domain);
        let err4 = solve(1000.0, &domain);

        assert!(err2 < err1);
        assert!(err3 < err2);
        assert!(err4 < err3);
    }


    #[test]
    fn sample_4() {
        /* square into triangles */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let e1 = Rc::new(Edge::new(&p1, &p2));
        let e2 = Rc::new(Edge::new(&p2, &p3));
        let e3 = Rc::new(Edge::new(&p3, &p4));
        let e4 = Rc::new(Edge::new(&p4, &p1));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p4, &p2, &p3));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e3, vec![1.0, 1.0]);
        domain.insert_neumann_constraint(&e2, vec![0.0, 0.0]);
        domain.insert_neumann_constraint(&e4, vec![0.0, 0.0]);

        fn solve(sigma: f64, domain: &Domain) -> f64 {
            let (system_matrix, extern_matrix) = build(sigma, &domain);
            let inverse = system_matrix.try_inverse().unwrap();
            let answer = inverse * extern_matrix;

            let expected =
                DMatrix::<f64>::from_row_slice(6, 1, &vec![0.0, 0.0, 1.0, 1.0, 0.0, 1.0]);
            let error = ((&answer - &expected).transpose() * (&answer - &expected))[(0, 0)].sqrt();
            println!("{} {}", answer, error);
            return error;
        }
        let err1 = solve(10.0, &domain);
        let err2 = solve(100.0, &domain);
        let err3 = solve(1000.0, &domain);
        let err4 = solve(10000.0, &domain);
        let err5 = solve(100000.0, &domain);

        assert!(err2 < err1);
        assert!(err3 < err2);
        assert!(err4 < err3);
        assert!(err5 < err4);
    }
}
