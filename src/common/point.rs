use float_cmp;
use nalgebra::Matrix2x1;
use num::Float;
use std::hash::{Hash, Hasher};

use std::fmt;

pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x: x, y: y }
    }

    pub fn as_matrix(&self) -> Matrix2x1<f64> {
        Matrix2x1::new(self.x, self.y)
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

/* Hash implementation */
impl Hash for Point {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let (m, e, s) = Float::integer_decode(self.x);
        m.hash(state);
        e.hash(state);
        s.hash(state);

        let (m, e, s) = Float::integer_decode(self.y);
        m.hash(state);
        e.hash(state);
        s.hash(state);
    }
}

/* Equality implementation */
impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        return float_cmp::approx_eq!(f64, self.x, other.x, epsilon = 1.0E-14f64)
            && float_cmp::approx_eq!(f64, self.y, other.y, epsilon = 1.0E-14f64);
    }
}

impl Eq for Point {}
