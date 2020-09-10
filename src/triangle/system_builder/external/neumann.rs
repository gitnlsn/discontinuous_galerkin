use nalgebra::DMatrix;

use crate::triangle::{integrands::neumann_constraint, system_builder::domain::Domain};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(b_matrix: &mut DMatrix<f64>, domain: &Domain) -> Result<(), ()> {
    let interfaces = domain.interfaces();
    for d_constraint in domain.neumann_constraints.iter() {
        let element = Rc::clone(&d_constraint.element);
        let inner_edge = Rc::clone(&d_constraint.boundary_edge);

        let p1 = Rc::clone(&inner_edge.p1);
        let p2 = Rc::clone(&inner_edge.p2);
        let p3 = Rc::clone(&element.opposite_vertex(&inner_edge).unwrap());

        let u1 = d_constraint.values.get(&p1).unwrap();
        let u2 = d_constraint.values.get(&p2).unwrap();

        let global_p1 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&p1)))
            .unwrap();

        let global_p2 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&p2)))
            .unwrap();

        let global_p3 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&p3)))
            .unwrap();

        let local_p1: usize = 0;
        let local_p2: usize = 1;
        let local_p3: usize = 2;

        let flux_matrix = neumann_constraint::neumann(&p1, &p2, &p3, *u1, *u2);

        /* Populating system_matrix upper left */
        b_matrix[(global_p1, 0)] += flux_matrix[(local_p1, 0)];
        b_matrix[(global_p2, 0)] += flux_matrix[(local_p2, 0)];
        b_matrix[(global_p3, 0)] += flux_matrix[(local_p3, 0)];
    }
    return Ok(());
}

#[cfg(test)]
mod build {
    use super::*;
    use crate::common::{point::Point, edge::Edge};
    use crate::triangle::element::TriangleElementL1;

    #[test]
    fn sample_1() {
        /* square into triangles */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let e1 = Rc::new(Edge::new(&p1, &p2));
        let e2 = Rc::new(Edge::new(&p3, &p4));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p4, &p2, &p3));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        domain.insert_neumann_constraint(&e1, vec![0.0, 0.0]);
        domain.insert_neumann_constraint(&e2, vec![2.0, 2.0]);

        let mut b_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(6, 1);

        build(&mut b_matrix, &domain);

        println!("{}", b_matrix);
    }
}
