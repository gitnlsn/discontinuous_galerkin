use crate::common::point::Point;
use crate::triangle::element::TriangleElementL1;

use std::collections::HashMap;
use std::rc::Rc;

/**
 * Abstraction of external fields acting on the project domain
 */
pub struct External {
    pub element: Rc<TriangleElementL1>,
    pub value: HashMap<Rc<Point>, f64>,
}

/**
 * Abstraction of variable fields representing the project domain
 */
pub struct Internal {
    pub element: Rc<TriangleElementL1>,
    pub value: HashMap<Rc<Point>, f64>,
}
