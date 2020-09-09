use crate::common::{edge::Edge, point::Point};

use crate::triangle::element::TriangleElementL1;

use std::collections::HashMap;
use std::rc::Rc;

pub struct BoundaryConstraint {
    pub element: Rc<TriangleElementL1>,
    pub boundary_edge: Rc<Edge>,
    pub values: HashMap<Rc<Point>, f64>,
}
