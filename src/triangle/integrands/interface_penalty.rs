use crate::common::point::Point;
use nalgebra::Matrix3;

use crate::triangle::integrands::utils;

/**
 * For the default L1 triangle, flux matrix only depends
 * on the points located at the interface
 */
pub fn bilinear_penalty(
    p1: &Point,
    p2: &Point,
    p3: &Point,
    p4: &Point,
    p5: &Point,
    p6: &Point,
    edge_index: usize,
) -> Matrix3<f64> {
    let sqrt_2 = (2.0 as f64).sqrt();
    let base = match edge_index {
        /* 0-edge (0,0) -> (1,0) */
        0 => Matrix3::new(
            -(p1.x + p2.x + p1.y + p2.y) / 2.0 + 1.0,
            -(p1.x / 6.0 + p2.x / 3.0 + p1.y / 6.0 + p2.y / 3.0) + 0.5,
            0.0,
            p1.x / 2.0 + p2.x / 2.0,
            p1.x / 6.0 + p2.x / 3.0,
            0.0,
            p1.y / 2.0 + p2.y / 2.0,
            p1.y / 6.0 + p2.y / 3.0,
            0.0,
        ),
        /* 1-edge (1,0) -> (0,1) */
        1 => Matrix3::new(
            (2.0 - p2.x - p3.x - p2.y - p3.y) * sqrt_2 / 2.0,
            (3.0 - 2.0 * p2.x - p3.x - 2.0 * p2.y - p3.y) * sqrt_2 / 6.0,
            (3.0 - p2.x - 2.0 * p3.x - p2.y - 2.0 * p3.y) * sqrt_2 / 6.0,
            (p2.x + p3.x) * sqrt_2 / 2.0,
            (2.0 * p2.x + p3.x) * sqrt_2 / 6.0,
            (p2.x + 2.0 * p3.x) * sqrt_2 / 6.0,
            (p2.y + p3.y) * sqrt_2 / 2.0,
            (2.0 * p2.y + p3.y) * sqrt_2 / 6.0,
            (p2.y + 2.0 * p3.y) * sqrt_2 / 6.0,
        ),
        /* 2-edge (0,1) -> (0,0) */
        2 => Matrix3::new(
            (p1.x + p3.x + p1.y + p3.y) / 2.0 - 1.0,
            0.0,
            p1.x / 6.0 + p3.x / 3.0 + p1.y / 6.0 + p3.y / 3.0 - 0.5,
            -p1.x / 2.0 - p3.x / 2.0,
            0.0,
            -p1.x / 6.0 - p3.x / 3.0,
            -p1.y / 2.0 - p3.y / 2.0,
            0.0,
            -p1.y / 6.0 - p3.y / 3.0,
        ),
        _ => panic!("Not expecting greater edge indices"),
    };
    return base
        * utils::coordinante_transformation(&p1, &p2, &p3).transpose()
        * utils::field_transformation(&p4, &p5, &p6)
            .try_inverse()
            .unwrap();
}

#[cfg(test)]
mod linear_natural {
    use super::*;
    use crate::triangle::integrands::dirichlet_constraint;
    use crate::triangle::system_builder::assembler_utils;
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
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        /* Assembling system */
        let mut system_matrix = DMatrix::<f64>::zeros(6, 6);
        let map_11 = assembler_utils::square_map(0, 1, 2);
        let map_12 = assembler_utils::cross_map(0, 1, 2, 3, 4, 5);
        let map_21 = assembler_utils::cross_map(3, 4, 5, 0, 1, 2);
        let map_22 = assembler_utils::square_map(3, 4, 5);

        assembler_utils::map(
            /* mapping (0,0)-(1,0) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p1, &p2, &p4, 0)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_11,
        );
        assembler_utils::map(
            /* mapping (0,1)-(0,0) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p1, &p2, &p4, 2)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_11,
        );
        assembler_utils::map(
            /* mapping (1,1)-(0,1) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p3, &p4, &p2, 0)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_22,
        );
        assembler_utils::map(
            /* mapping (1,0)-(1,1) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p3, &p4, &p2, 2)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_22,
        );
        /* Interface terms */
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut system_matrix,
            &bilinear_penalty(&p1, &p2, &p4, &p1, &p2, &p4, 1)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_11,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut system_matrix,
            &(-bilinear_penalty(&p1, &p2, &p4, &p3, &p4, &p2, 1))
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_12,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut system_matrix,
            &bilinear_penalty(&p3, &p4, &p2, &p3, &p4, &p2, 1)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_22,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut system_matrix,
            &(-bilinear_penalty(&p3, &p4, &p2, &p1, &p2, &p4, 1))
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &map_21,
        );

        // println!("{}", system_matrix);

        /* Extern matrix */
        let mut extern_matrix = DMatrix::<f64>::zeros(6, 1);
        let map_1 = assembler_utils::linear_map(0, 1, 2);
        let map_2 = assembler_utils::linear_map(3, 4, 5);
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p1, &p2, &p4, 0.0, 0.0, 2.0, 0)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &map_1,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p1, &p2, &p4, 0.0, 0.0, 2.0, 2)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &map_1,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p3, &p4, &p2, 0.0, 2.0, 0.0, 0)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &map_2,
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p3, &p4, &p2, 0.0, 2.0, 0.0, 2)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &map_2,
        );

        let solution = system_matrix.try_inverse().unwrap() * extern_matrix;
        let expected = DMatrix::<f64>::from_row_slice(6, 1, &vec![0.0, 0.0, 2.0, 0.0, 2.0, 0.0]);

        let error = ((&solution - &expected).transpose() * (&solution - &expected))[(0,0)].sqrt();
        // println!("{}", solution);
        // println!("{}", expected);
        // println!("{}", error);

        assert!(error < 1.0);/* Do not trust this test */
    }
}
