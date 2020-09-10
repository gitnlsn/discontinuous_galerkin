use nalgebra::DMatrix;

use crate::triangle::{integrands::penalty_zero, system_builder::domain::Domain};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(system_matrix: &mut DMatrix<f64>, sigma: f64, domain: &Domain) -> Result<(), ()> {
    let interfaces = domain.interfaces();
    for (e1, e2) in interfaces.iter() {
        let t_left = Rc::clone(domain.adjacency.get(e1).unwrap());
        let t_right = Rc::clone(domain.adjacency.get(e2).unwrap());

        let e1_opposed = t_left.opposite_vertex(e1).unwrap();
        let e2_opposed = t_right.opposite_vertex(e2).unwrap();

        let flux_matrix_left = penalty_zero::half_penalty(&e1.p1, &e1.p2) * sigma;
        let flux_matrix_right = penalty_zero::half_penalty(&e2.p1, &e2.p2) * sigma;

        let left_p1 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_left), Rc::clone(&e1.p1)))
            .unwrap();
        let left_p2 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_left), Rc::clone(&e1.p2)))
            .unwrap();
        let left_p3 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_left), Rc::clone(&e1_opposed)))
            .unwrap();

        let right_p1 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_right), Rc::clone(&e2.p1)))
            .unwrap();
        let right_p2 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_right), Rc::clone(&e2.p2)))
            .unwrap();
        let right_p3 = *domain
            .index_mapping
            .get(&(Rc::clone(&t_right), Rc::clone(&e2_opposed)))
            .unwrap();

        let lo_p1: usize = 0;
        let lo_p2: usize = 1;
        let lo_p3: usize = 2;

        /* Populating system_matrix upper left */
        system_matrix[(left_p1, left_p1)] += flux_matrix_left[(lo_p1, lo_p1)];
        system_matrix[(left_p1, left_p2)] += flux_matrix_left[(lo_p1, lo_p2)];
        system_matrix[(left_p1, left_p3)] += flux_matrix_left[(lo_p1, lo_p3)];

        system_matrix[(left_p2, left_p1)] += flux_matrix_left[(lo_p2, lo_p1)];
        system_matrix[(left_p2, left_p2)] += flux_matrix_left[(lo_p2, lo_p2)];
        system_matrix[(left_p2, left_p3)] += flux_matrix_left[(lo_p2, lo_p3)];

        system_matrix[(left_p3, left_p1)] += flux_matrix_left[(lo_p3, lo_p1)];
        system_matrix[(left_p3, left_p2)] += flux_matrix_left[(lo_p3, lo_p2)];
        system_matrix[(left_p3, left_p3)] += flux_matrix_left[(lo_p3, lo_p3)];

        /* Populating system_matrix upper right */
        system_matrix[(left_p1, right_p1)] -= flux_matrix_left[(lo_p1, lo_p1)];
        system_matrix[(left_p1, right_p2)] -= flux_matrix_left[(lo_p1, lo_p2)];
        system_matrix[(left_p1, right_p3)] -= flux_matrix_left[(lo_p1, lo_p3)];

        system_matrix[(left_p2, right_p1)] -= flux_matrix_left[(lo_p2, lo_p1)];
        system_matrix[(left_p2, right_p2)] -= flux_matrix_left[(lo_p2, lo_p2)];
        system_matrix[(left_p2, right_p3)] -= flux_matrix_left[(lo_p2, lo_p3)];

        system_matrix[(left_p3, right_p1)] -= flux_matrix_left[(lo_p3, lo_p1)];
        system_matrix[(left_p3, right_p2)] -= flux_matrix_left[(lo_p3, lo_p2)];
        system_matrix[(left_p3, right_p3)] -= flux_matrix_left[(lo_p3, lo_p3)];

        /* Populating system_matrix lower left */
        system_matrix[(right_p1, left_p1)] -= flux_matrix_right[(lo_p1, lo_p1)];
        system_matrix[(right_p1, left_p2)] -= flux_matrix_right[(lo_p1, lo_p2)];
        system_matrix[(right_p1, left_p3)] -= flux_matrix_right[(lo_p1, lo_p3)];

        system_matrix[(right_p2, left_p1)] -= flux_matrix_right[(lo_p2, lo_p1)];
        system_matrix[(right_p2, left_p2)] -= flux_matrix_right[(lo_p2, lo_p2)];
        system_matrix[(right_p2, left_p3)] -= flux_matrix_right[(lo_p2, lo_p3)];

        system_matrix[(right_p3, left_p1)] -= flux_matrix_right[(lo_p3, lo_p1)];
        system_matrix[(right_p3, left_p2)] -= flux_matrix_right[(lo_p3, lo_p2)];
        system_matrix[(right_p3, left_p3)] -= flux_matrix_right[(lo_p3, lo_p3)];

        /* Populating system_matrix lower right */
        system_matrix[(right_p1, right_p1)] += flux_matrix_right[(lo_p1, lo_p1)];
        system_matrix[(right_p1, right_p2)] += flux_matrix_right[(lo_p1, lo_p2)];
        system_matrix[(right_p1, right_p3)] += flux_matrix_right[(lo_p1, lo_p3)];

        system_matrix[(right_p2, right_p1)] += flux_matrix_right[(lo_p2, lo_p1)];
        system_matrix[(right_p2, right_p2)] += flux_matrix_right[(lo_p2, lo_p2)];
        system_matrix[(right_p2, right_p3)] += flux_matrix_right[(lo_p2, lo_p3)];

        system_matrix[(right_p3, right_p1)] += flux_matrix_right[(lo_p3, lo_p1)];
        system_matrix[(right_p3, right_p2)] += flux_matrix_right[(lo_p3, lo_p2)];
        system_matrix[(right_p3, right_p3)] += flux_matrix_right[(lo_p3, lo_p3)];
    }
    return Ok(());
}

#[cfg(test)]
mod build {
    use super::*;
    use crate::common::point::Point;
    use crate::triangle::element::TriangleElementL1;

    #[test]
    fn sample_1() {
        /* square into triangles */
        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(1.0, 1.0));
        let p4 = Rc::new(Point::new(0.0, 1.0));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p4, &p2, &p3));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        let mut system_matrix: DMatrix<f64> = DMatrix::<f64>::zeros(6, 6);

        build(&mut system_matrix, 3.0, &domain);

        println!("{}", system_matrix);
    }

    #[test]
    fn sample_2() {
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

        let variables_length = domain.elements.len() * 3;
        let mut system_matrix: DMatrix<f64> =
            DMatrix::<f64>::zeros(variables_length, variables_length);

        build(&mut system_matrix, 3.0, &domain);

        println!("{}", system_matrix);
    }
}
