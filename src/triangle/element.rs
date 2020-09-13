use crate::common::edge::Edge;
use crate::common::point::Point;

use std::fmt;
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

impl fmt::Display for TriangleElementL1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {}", self.p1, self.p2, self.p3)
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

    pub fn edge_index(&self, edge: &Edge) -> Result<usize, ()> {
        if edge.p1 == self.p1 && edge.p2 == self.p2 {
            return Ok(0);
        } else if edge.p1 == self.p2 && edge.p2 == self.p3 {
            return Ok(1);
        } else if edge.p1 == self.p3 && edge.p2 == self.p1 {
            return Ok(2);
        }
        return Err(());
    }
}

#[cfg(test)]
mod hash {
    use super::*;

    #[test]
    fn sample_1() {
        use std::collections::HashSet;

        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(0.0, 1.0));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p3));
        let t2 = Rc::new(TriangleElementL1::new(&p2, &p3, &p1));

        let mut set: HashSet<Rc<TriangleElementL1>> = HashSet::new();

        set.insert(Rc::clone(&t1));
        set.insert(Rc::clone(&t2));
        assert_eq!(set.len(), 2);

        set.insert(Rc::clone(&t1));
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn sample_2() {
        use std::collections::HashSet;

        let p1 = Rc::new(Point::new(0.0, 0.0));
        let p2 = Rc::new(Point::new(1.0, 0.0));
        let p3 = Rc::new(Point::new(0.0, 1.0));

        let t1 = Rc::new(TriangleElementL1::new(&p1, &p2, &p3));
        let t2 = Rc::new(TriangleElementL1::new(&p2, &p3, &p1));

        let mut set: HashSet<(Rc<TriangleElementL1>, Rc<Point>)> = HashSet::new();

        set.insert((Rc::clone(&t1), Rc::clone(&p1)));
        set.insert((Rc::clone(&t2), Rc::clone(&p1)));
        assert_eq!(set.len(), 2);

        set.insert((Rc::clone(&t1), Rc::clone(&p1)));
        assert_eq!(set.len(), 2);
    }
}
