use nalgebra::DMatrix;

use crate::triangle::{integrands::mass, system_builder::domain::Domain};

use std::rc::Rc;

/**
 * Fills the system matrix with mass matrix according to each element
 */
pub fn build(system_matrix: &mut DMatrix<f64>, domain: &Domain) -> Result<(), ()> {
    for index in 0..domain.elements.len() {
        let element = Rc::clone(domain.elements.get(index).unwrap());

        let gl_p1 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p1)))
            .unwrap();
        let gl_p2 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p2)))
            .unwrap();
        let gl_p3 = *domain
            .index_mapping
            .get(&(Rc::clone(&element), Rc::clone(&element.p3)))
            .unwrap();

        let mass_matrix = mass::matrix(&element.p1, &element.p2, &element.p3);

        let lo_p1: usize = 0;
        let lo_p2: usize = 1;
        let lo_p3: usize = 2;

        system_matrix[(gl_p1, gl_p1)] += mass_matrix[(lo_p1, lo_p1)];
        system_matrix[(gl_p1, gl_p2)] += mass_matrix[(lo_p1, lo_p2)];
        system_matrix[(gl_p1, gl_p3)] += mass_matrix[(lo_p1, lo_p3)];

        system_matrix[(gl_p2, gl_p1)] += mass_matrix[(lo_p2, lo_p1)];
        system_matrix[(gl_p2, gl_p2)] += mass_matrix[(lo_p2, lo_p2)];
        system_matrix[(gl_p2, gl_p3)] += mass_matrix[(lo_p2, lo_p3)];

        system_matrix[(gl_p3, gl_p1)] += mass_matrix[(lo_p3, lo_p1)];
        system_matrix[(gl_p3, gl_p2)] += mass_matrix[(lo_p3, lo_p2)];
        system_matrix[(gl_p3, gl_p3)] += mass_matrix[(lo_p3, lo_p3)];
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

        build(&mut system_matrix, &domain);

        assert_eq!(system_matrix[(0, 0)], 2.0);
        assert_eq!(system_matrix[(0, 1)], -1.0);
        assert_eq!(system_matrix[(0, 2)], -1.0);
        assert_eq!(system_matrix[(1, 1)], 1.0);
        assert_eq!(system_matrix[(2, 2)], 1.0);
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

        build(&mut system_matrix, &domain);

        assert_eq!(system_matrix[(0, 0)], 2.0);
        assert_eq!(system_matrix[(0, 1)], -2.0);
        assert_eq!(system_matrix[(0, 2)], 0.0);
        assert_eq!(system_matrix[(1, 1)], 2.5);
        assert_eq!(system_matrix[(1, 2)], -0.5);
        assert_eq!(system_matrix[(2, 2)], 0.5);
    }
}
