use crate::common::point::Point;
use std::hash::Hash;
use std::rc::Rc;
use std::collections::HashSet;

#[derive(Hash)]
pub struct Edge {
    pub p1: Rc<Point>,
    pub p2: Rc<Point>,
}

impl Edge {
    pub fn new(p1: &Rc<Point>, p2: &Rc<Point>) -> Self {
        Self {
            p1: Rc::clone(p1),
            p2: Rc::clone(p2),
        }
    }

    pub fn opposed(&self) -> Self {
        Self {
            p1: Rc::clone(&self.p2),
            p2: Rc::clone(&self.p1),
        }
    }

    pub fn double_edges(edges: Vec<Rc<Edge>>) -> Vec<(Rc<Edge>, Rc<Edge>)> {
        let mut aux_list: HashSet<Rc<Edge>> = HashSet::new();
        let mut double_edge_list: Vec<(Rc<Edge>, Rc<Edge>)> = Vec::new();

        for edge in edges.iter() {
            if aux_list.contains(&edge.opposed()) {
                aux_list.remove(edge);
                let edge = Rc::clone(edge);
                let opposed = Rc::new(edge.opposed());
                double_edge_list.push((edge, opposed));
            } else {
                let edge = Rc::clone(edge);
                aux_list.insert(edge);
            }
        }

        return double_edge_list;
    }
}

/* Equality implementation */
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        return self.p1 == other.p1 && self.p2 == other.p2;
    }
}

impl Eq for Edge {}
