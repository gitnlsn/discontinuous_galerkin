use nalgebra::DMatrix;

use crate::triangle::{integrands::dirichlet_constraint, system_builder::domain::Domain};

use std::rc::Rc;

/**
 *  Fills both matrices
 */
pub fn build(
    system_matrix: &mut DMatrix<f64>, /* NxN matrix */
    extern_matrix: &mut DMatrix<f64>, /* Nx1 matrix */
    domain: &Domain,
) -> Result<(), ()> {
    for d_constraint in domain.dirichlet_constraints.iter() {
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

        /* Bilinear Penalty */
        let penalty_bilinear = dirichlet_constraint::dirichlet_bilinear_penalty(&p1, &p2, &p3);
        system_matrix[(global_p1, global_p1)] += penalty_bilinear[(local_p1, local_p1)];
        system_matrix[(global_p1, global_p2)] += penalty_bilinear[(local_p1, local_p2)];
        system_matrix[(global_p1, global_p3)] += penalty_bilinear[(local_p1, local_p3)];

        system_matrix[(global_p2, global_p1)] += penalty_bilinear[(local_p2, local_p1)];
        system_matrix[(global_p2, global_p2)] += penalty_bilinear[(local_p2, local_p2)];
        system_matrix[(global_p2, global_p3)] += penalty_bilinear[(local_p2, local_p3)];

        system_matrix[(global_p3, global_p1)] += penalty_bilinear[(local_p3, local_p1)];
        system_matrix[(global_p3, global_p2)] += penalty_bilinear[(local_p3, local_p2)];
        system_matrix[(global_p3, global_p3)] += penalty_bilinear[(local_p3, local_p3)];

        /* Linear Natural  */
        let b1_matrix =
            dirichlet_constraint::dirichlet_linear_natural(&p1, &p2, &p3, *u1, *u2);

        extern_matrix[(global_p1, 0)] += b1_matrix[(local_p1, 0)];
        extern_matrix[(global_p2, 0)] += b1_matrix[(local_p2, 0)];
        extern_matrix[(global_p3, 0)] += b1_matrix[(local_p3, 0)];

        /* Linear Penalty */
        let penalty_linear =
            dirichlet_constraint::dirichlet_linear_penalty(&p1, &p2, &p3, *u1, *u2);

        extern_matrix[(global_p1, 0)] += penalty_linear[(local_p1, 0)];
        extern_matrix[(global_p2, 0)] += penalty_linear[(local_p2, 0)];
        extern_matrix[(global_p3, 0)] += penalty_linear[(local_p3, 0)];
    }

    return Ok(());
}

#[cfg(test)]
mod dirichlet {
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

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e2, vec![2.0, 2.0]);

        let mut system_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(6, 6);
        let mut extern_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(6, 1);

        build(&mut system_matrix, &mut extern_matrix, &domain);

        println!("{}", system_matrix);
        println!("{}", extern_matrix);
    }
}