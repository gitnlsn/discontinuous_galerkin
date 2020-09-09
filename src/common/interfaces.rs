use crate::common::point::Point;
use crate::triangle::element::TriangleElementL1;

use std::rc::Rc;

pub struct ElementsInterface {
    pub e1: Rc<TriangleElementL1>,
    pub e2: Rc<TriangleElementL1>,
}

pub struct OuterInterface {
    pub e1: Rc<TriangleElementL1>,
    pub edge: (Rc<Point>, Rc<Point>),
}
