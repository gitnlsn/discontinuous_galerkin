use nalgebra::DMatrix;

use crate::triangle::{
    integrands::flux_artificial,
    system_builder::{assembler_utils, domain::Domain},
};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(system_matrix: &mut DMatrix<f64>, domain: &Domain) -> Result<(), ()> {
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
                    &flux_artificial::half_flux(
                        &t_left.p1,
                        &t_left.p2,
                        &t_left.p3,
                        &t_left.p1,
                        &t_left.p2,
                        &t_left.p3,
                        t_left.edge_index(&edge).unwrap(),
                    )
                    .slice((0, 0), (3, 3))
                    .clone_owned(),
                    &assembler_utils::square_map(left_p1, left_p2, left_p3),
                );

                assembler_utils::map(
                    /* mapping inner interface (0,1)-(1,0) */
                    system_matrix,
                    &flux_artificial::half_flux(
                        &t_left.p1,
                        &t_left.p2,
                        &t_left.p3,
                        &t_right.p1,
                        &t_right.p2,
                        &t_right.p3,
                        t_left.edge_index(&edge).unwrap(),
                    )
                    .slice((0, 0), (3, 3))
                    .clone_owned(),
                    &assembler_utils::cross_map(
                        left_p1, left_p2, left_p3, right_p1, right_p2, right_p3,
                    ),
                );
            } /* end - adjacency */
        } /* end - for edge in triangle */
    } /* end - for element in domain */
    return Ok(());
}
