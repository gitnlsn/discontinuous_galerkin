use crate::common::edge::Edge;
use crate::common::point::Point;
use crate::triangle::{boundary_constraint::BoundaryConstraint, element::TriangleElementL1};

use std::collections::HashMap;
use std::rc::Rc;

pub struct Domain {
    pub elements: Vec<Rc<TriangleElementL1>>,
    pub adjacency: HashMap<Rc<Edge>, Rc<TriangleElementL1>>,

    /*
        Index mapping will point each triangle point directly into a index
        Indexation is done sequentially, starting from 0
    */
    pub index_mapping: HashMap<(Rc<TriangleElementL1>, Rc<Point>), usize>,

    /*
        Boundary Constraints
    */
    pub dirichlet_constraints: Vec<Rc<BoundaryConstraint>>,
    pub neumann_constraints: Vec<Rc<BoundaryConstraint>>,
}

impl Domain {
    pub fn new_empty() -> Self {
        Domain {
            elements: Vec::new(),
            adjacency: HashMap::new(),
            index_mapping: HashMap::new(),
            dirichlet_constraints: Vec::new(),
            neumann_constraints: Vec::new(),
        }
    }

    pub fn insert_element(&mut self, triangle: &Rc<TriangleElementL1>) {
        let (e1, e2, e3) = triangle.inner_edges();
        self.elements.push(Rc::clone(triangle));

        self.adjacency.insert(Rc::clone(&e1), Rc::clone(triangle));
        self.adjacency.insert(Rc::clone(&e2), Rc::clone(triangle));
        self.adjacency.insert(Rc::clone(&e3), Rc::clone(triangle));

        self.index_mapping.insert(
            (Rc::clone(triangle), Rc::clone(&triangle.p1)),
            self.index_mapping.len(),
        );
        self.index_mapping.insert(
            (Rc::clone(triangle), Rc::clone(&triangle.p2)),
            self.index_mapping.len(),
        );
        self.index_mapping.insert(
            (Rc::clone(triangle), Rc::clone(&triangle.p3)),
            self.index_mapping.len(),
        );
    }

    pub fn insert_dirichlet_constraint(&mut self, edge: &Rc<Edge>, values: Vec<f64>) {
        let triangle = self.adjacency.get(edge).unwrap();
        let mut values_mapping: HashMap<Rc<Point>, f64> = HashMap::new();

        values_mapping.insert(Rc::clone(&edge.p1), *values.get(0).unwrap());
        values_mapping.insert(Rc::clone(&edge.p2), *values.get(1).unwrap());

        self.dirichlet_constraints.push(Rc::new(BoundaryConstraint {
            element: Rc::clone(triangle),
            boundary_edge: Rc::clone(edge),
            values: values_mapping,
        }))
    }

    pub fn insert_neumann_constraint(&mut self, edge: &Rc<Edge>, values: Vec<f64>) {
        let triangle = self.adjacency.get(edge).unwrap();
        let mut values_mapping: HashMap<Rc<Point>, f64> = HashMap::new();

        values_mapping.insert(Rc::clone(&edge.p1), *values.get(0).unwrap());
        values_mapping.insert(Rc::clone(&edge.p2), *values.get(1).unwrap());

        self.neumann_constraints.push(Rc::new(BoundaryConstraint {
            element: Rc::clone(triangle),
            boundary_edge: Rc::clone(edge),
            values: values_mapping,
        }))
    }
} /* end - domain */
