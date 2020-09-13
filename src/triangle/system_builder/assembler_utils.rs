use nalgebra::DMatrix;

use std::collections::HashMap;

fn linear(global_1: usize, global_2: usize, global_3: usize) -> HashMap<usize, usize> {
    let mut map: HashMap<usize, usize> = HashMap::new();

    map.insert(0, global_1);
    map.insert(1, global_2);
    map.insert(2, global_3);

    return map;
}

pub fn linear_map(
    global_1: usize,
    global_2: usize,
    global_3: usize,
) -> HashMap<(usize, usize), (usize, usize)> {
    let linear_map = linear(global_1, global_2, global_3);
    let mut squared_map: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

    for (local_row, global_row) in linear_map.iter() {
        squared_map.insert((*local_row, 0), (*global_row, 0));
    }

    return squared_map;
}

pub fn cross_map(
    global_1: usize,
    global_2: usize,
    global_3: usize,
    global_4: usize,
    global_5: usize,
    global_6: usize,
) -> HashMap<(usize, usize), (usize, usize)> {
    let l1 = linear(global_1, global_2, global_3);
    let l2 = linear(global_4, global_5, global_6);
    let mut squared_map: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

    for (local_row, global_row) in l1.iter() {
        for (local_col, global_col) in l2.iter() {
            squared_map.insert((*local_row, *local_col), (*global_row, *global_col));
        }
    }

    return squared_map;
}

pub fn square_map(
    global_1: usize,
    global_2: usize,
    global_3: usize,
) -> HashMap<(usize, usize), (usize, usize)> {
    let linear_map = linear(global_1, global_2, global_3);
    let mut squared_map: HashMap<(usize, usize), (usize, usize)> = HashMap::new();

    for (local_row, global_row) in linear_map.iter() {
        for (local_col, global_col) in linear_map.iter() {
            squared_map.insert((*local_row, *local_col), (*global_row, *global_col));
        }
    }

    return squared_map;
}

pub fn map(
    global: &mut DMatrix<f64>,
    local: &DMatrix<f64>,
    index_map: &HashMap<(usize, usize), (usize, usize)>,
) {
    for (local_index, global_index) in index_map.iter() {
        global[*global_index] += local[*local_index];
    }
}
