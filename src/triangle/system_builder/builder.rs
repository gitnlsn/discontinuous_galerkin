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

    mass::build(&mut system_matrix, domain);
    flux_natural::build(&mut system_matrix, domain);
    flux_artificial::build(&mut system_matrix, domain);

    jump_penalty::build(&mut system_matrix, sigma, domain);
    dirichlet::build(&mut system_matrix, &mut extern_matrix, sigma, domain);
    neumann::build(&mut extern_matrix, domain);

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

        let (system_matrix, extern_matrix) = build(10000.0, &domain);
        println!("{}", system_matrix);
        // println!("{}", extern_matrix);
        
        let inverse = system_matrix.try_inverse().unwrap();
        // println!("{}", inverse);

        let answer = inverse * extern_matrix;
        println!("{}", answer);
    }

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
        
        // let ed = Rc::new(Edge::new(&p2, &p4));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p4));
        let t2 = Rc::new(TriangleElementL1::new(&p4, &p2, &p3));

        let mut domain = Domain::new_empty();
        domain.insert_element(&t1);
        domain.insert_element(&t2);

        domain.insert_dirichlet_constraint(&e1, vec![0.0, 1.0]);
        // domain.insert_dirichlet_constraint(&e4, vec![0.0, 1.0]);
        // domain.insert_dirichlet_constraint(&ed, vec![0.0, 0.0]);
        
        domain.insert_dirichlet_constraint(&e2, vec![1.0, 1.0]);
        domain.insert_dirichlet_constraint(&e3, vec![0.0, 0.0]);
        domain.insert_dirichlet_constraint(&e4, vec![0.0, 0.0]);

        // domain.insert_neumann_constraint(&e2, vec![0.0, 0.0]);
        // domain.insert_neumann_constraint(&ed, vec![0.0, 0.0]);
        // domain.insert_neumann_constraint(&e4, vec![0.0, 0.0]);

        let (system_matrix, extern_matrix) = build(10000.0, &domain);
        println!("{}", system_matrix);
        // println!("{}", extern_matrix);

        let inverse = system_matrix.try_inverse().unwrap();
        // println!("{}", inverse);

        let answer = inverse * extern_matrix;
        println!("{}", answer);
    }
}
