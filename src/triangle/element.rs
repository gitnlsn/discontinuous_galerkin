use crate::common::edge::Edge;
use crate::common::point::Point;

use std::rc::Rc;

#[derive(Hash)]
pub struct TriangleElementL1 {
    pub p1: Rc<Point>,
    pub p2: Rc<Point>,
    pub p3: Rc<Point>,
}

/* Equality implementation */
impl PartialEq for TriangleElementL1 {
    fn eq(&self, other: &Self) -> bool {
        return self.p1 == other.p1 && self.p2 == other.p2 && self.p3 == other.p3;
    }
}

impl TriangleElementL1 {
    pub fn opposite_vertex(&self, edge: &Edge) -> Result<Rc<Point>, ()> {
        if self.p1 == edge.p1 {
            return Ok(Rc::clone(&self.p3));
        } else if self.p2 == edge.p1 {
            return Ok(Rc::clone(&self.p1));
        } else if self.p3 == edge.p1 {
            return Ok(Rc::clone(&self.p2));
        } else {
            return Err(());
        }
    }
}

impl Eq for TriangleElementL1 {}

impl TriangleElementL1 {
    pub fn new(p1: &Rc<Point>, p2: &Rc<Point>, p3: &Rc<Point>) -> Self {
        Self {
            p1: Rc::clone(p1),
            p2: Rc::clone(p2),
            p3: Rc::clone(p3),
        }
    }

    pub fn inner_edges(&self) -> (Rc<Edge>, Rc<Edge>, Rc<Edge>) {
        let e1 = Rc::new(Edge::new(&self.p1, &self.p2));
        let e2 = Rc::new(Edge::new(&self.p2, &self.p3));
        let e3 = Rc::new(Edge::new(&self.p3, &self.p1));
        return (e1, e2, e3);
    }
}
