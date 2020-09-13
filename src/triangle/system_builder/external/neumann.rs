use nalgebra::DMatrix;

use crate::triangle::{
    integrands::neumann_constraint,
    system_builder::{assembler_utils, domain::Domain},
};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(b_matrix: &mut DMatrix<f64>, domain: &Domain) -> Result<(), ()> {
    for d_constraint in domain.neumann_constraints.iter() {
        let element = Rc::clone(&d_constraint.element);
        let inner_edge = Rc::clone(&d_constraint.boundary_edge);

        let u1: f64 = match d_constraint.values.get(&element.p1) {
            Some(value) => *value,
            None => 0.0,
        };
        let u2: f64 = match d_constraint.values.get(&element.p2) {
            Some(value) => *value,
            None => 0.0,
        };
        let u3: f64 = match d_constraint.values.get(&element.p3) {
            Some(value) => *value,
            None => 0.0,
        };

        let global_p1 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p1)))
            .unwrap();

        let global_p2 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p2)))
            .unwrap();

        let global_p3 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p3)))
            .unwrap();

        assembler_utils::map(
            b_matrix,
            &neumann_constraint::neumann(
                &element.p1,
                &element.p2,
                &element.p3,
                u1,
                u2,
                u3,
                element.edge_index(&inner_edge).unwrap(),
            )
            .slice((0, 0), (3, 1))
            .clone_owned(),
            &assembler_utils::linear_map(global_p1, global_p2, global_p3),
        );
    }
    return Ok(());
}
