use nalgebra::DMatrix;

use crate::triangle::{
    integrands::interface_penalty,
    system_builder::{assembler_utils, domain::Domain},
};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(system_matrix: &mut DMatrix<f64>, sigma: f64, domain: &Domain) -> Result<(), ()> {
    for element in domain.elements.iter() {
        let (e1, e2, e3) = element.inner_edges();
        for edge in vec![e1, e2, e3].iter() {
            let edge = Rc::clone(edge);
            let edge_opp = Rc::new(edge.opposed());
            if domain.adjacency.contains_key(&edge_opp) {
                let t_left = Rc::clone(domain.adjacency.get(&edge).unwrap());
                let t_right = Rc::clone(domain.adjacency.get(&edge_opp).unwrap());

                let left_p1 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_left), Rc::clone(&t_left.p1)))
                    .unwrap();
                let left_p2 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_left), Rc::clone(&t_left.p2)))
                    .unwrap();
                let left_p3 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_left), Rc::clone(&t_left.p3)))
                    .unwrap();

                let right_p1 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_right), Rc::clone(&t_right.p1)))
                    .unwrap();
                let right_p2 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_right), Rc::clone(&t_right.p2)))
                    .unwrap();
                let right_p3 = *domain
                    .index_mapping
                    .get(&(Rc::clone(&t_right), Rc::clone(&t_right.p3)))
                    .unwrap();

                assembler_utils::map(
                    /* mapping inner interface (0,1)-(1,0) */
                    system_matrix,
                    &(sigma
                        * interface_penalty::bilinear_penalty(
                            &t_left.p1,
                            &t_left.p2,
                            &t_left.p3,
                            &t_left.p1,
                            &t_left.p2,
                            &t_left.p3,
                            t_left.edge_index(&edge).unwrap(),
                        ))
                    .slice((0, 0), (3, 3))
                    .clone_owned(),
                    &assembler_utils::square_map(left_p1, left_p2, left_p3),
                );

                assembler_utils::map(
                    /* mapping inner interface (0,1)-(1,0) */
                    system_matrix,
                    &(-sigma
                        * interface_penalty::bilinear_penalty(
                            &t_left.p1,
                            &t_left.p2,
                            &t_left.p3,
                            &t_right.p1,
                            &t_right.p2,
                            &t_right.p3,
                            t_left.edge_index(&edge).unwrap(),
                        ))
                    .slice((0, 0), (3, 3))
                    .clone_owned(),
                    &assembler_utils::cross_map(
                        left_p1, left_p2, left_p3, right_p1, right_p2, right_p3,
                    ),
                );
            } /* end - if adjancency */
        } /* end - for edge in triangle */
    } /* end - for element in domain */
    return Ok(());
}

#[cfg(test)]
mod build {
    use super::*;
    use crate::common::point::Point;
    use crate::triangle::{
        element::TriangleElementL1, integrands::dirichlet_constraint,
        system_builder::assembler_utils,
    };

    #[test]
    fn sample_1() {
        /*
            Weak imposed continuity is performed in triangle::integrand
            This test asserts system build method will result in the same error.
        */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p3, &p4, &p2));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        let mut system_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(6, 6);

        assembler_utils::map(
            /* mapping (0,0)-(1,0) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p1, &p2, &p4, 0)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &assembler_utils::square_map(0, 1, 2),
        );
        assembler_utils::map(
            /* mapping (0,1)-(0,0) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p1, &p2, &p4, 2)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &assembler_utils::square_map(0, 1, 2),
        );
        assembler_utils::map(
            /* mapping (1,1)-(0,1) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p3, &p4, &p2, 0)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &assembler_utils::square_map(3, 4, 5),
        );
        assembler_utils::map(
            /* mapping (1,0)-(1,1) */
            &mut system_matrix,
            &dirichlet_constraint::dirichlet_bilinear_penalty(&p3, &p4, &p2, 2)
                .slice((0, 0), (3, 3))
                .clone_owned(),
            &assembler_utils::square_map(3, 4, 5),
        );

        build(&mut system_matrix, 1.0, &domain);

        let mut extern_matrix = DMatrix::<f64>::zeros(6, 1);
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p1, &p2, &p4, 0.0, 0.0, 2.0, 0)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &assembler_utils::linear_map(0, 1, 2),
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p1, &p2, &p4, 0.0, 0.0, 2.0, 2)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &assembler_utils::linear_map(0, 1, 2),
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p3, &p4, &p2, 0.0, 2.0, 0.0, 0)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &assembler_utils::linear_map(3, 4, 5),
        );
        assembler_utils::map(
            /* mapping inner interface (0,1)-(1,0) */
            &mut extern_matrix,
            &dirichlet_constraint::dirichlet_linear_penalty(&p3, &p4, &p2, 0.0, 2.0, 0.0, 2)
                .slice((0, 0), (3, 1))
                .clone_owned(),
            &assembler_utils::linear_map(3, 4, 5),
        );

        let solution = system_matrix.try_inverse().unwrap() * extern_matrix;
        let expected = DMatrix::<f64>::from_row_slice(6, 1, &vec![0.0, 0.0, 2.0, 0.0, 2.0, 0.0]);

        let error = ((&solution - &expected).transpose() * (&solution - &expected))[(0, 0)].sqrt();

        assert_eq!(error, 0.375);
    }
}
