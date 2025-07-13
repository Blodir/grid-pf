use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
};

use ordered_float::OrderedFloat;

pub type Tile = (u32, u32);

// includes both start and destination elements
fn reconstruct_path(came_from: &HashMap<Tile, Tile>, mut current: Tile) -> Vec<Tile> {
    let mut total_path = vec![current];
    while let Some(&prev) = came_from.get(&current) {
        current = prev;
        total_path.push(current);
    }
    total_path.reverse();
    total_path
}

fn h(a: Tile, b: Tile) -> f32 {
    //(a.0.abs_diff(b.0) + a.1.abs_diff(b.1)) as f32
    ((b.0 as f32 - a.0 as f32).powi(2) + (b.1 as f32 - a.1 as f32).powi(2)).sqrt()
}

fn find_neighbors(
    current: &Tile,
    pathable: &Vec<Vec<bool>>,
    min_x: u32,
    max_x: u32,
    min_y: u32,
    max_y: u32,
) -> Vec<Tile> {
    let mut neighbors = vec![];
    let (x, y) = (current.0 as i32, current.1 as i32);

    let deltas = [
        (-1, 0),
        (1, 0),
        (0, -1),
        (0, 1),
        (-1, -1),
        (1, -1),
        (-1, 1),
        (1, 1),
    ];

    for (dx, dy) in deltas {
        let nx = x + dx;
        let ny = y + dy;

        if nx >= min_x as i32
            && nx <= max_x as i32
            && ny >= min_y as i32
            && ny <= max_y as i32
            && pathable[ny as usize][nx as usize]
        {
            neighbors.push((nx as u32, ny as u32));
        }
    }

    neighbors
}

pub fn astar(
    start: Tile,
    goal: Tile,
    pathable: &Vec<Vec<bool>>,
    min_x: u32,
    max_x: u32,
    min_y: u32,
    max_y: u32,
) -> Option<Vec<Tile>> {
    let mut came_from = HashMap::<Tile, Tile>::new();

    let mut g_score = HashMap::<Tile, f32>::new();
    g_score.insert(start, 0f32);

    let mut open_set = BinaryHeap::<(Reverse<OrderedFloat<f32>>, Tile)>::new();
    open_set.push((Reverse(OrderedFloat::from(h(start, goal))), start));

    while !open_set.is_empty() {
        let current = open_set.pop().unwrap().1;
        if current == goal {
            return Some(reconstruct_path(&came_from, current));
        }

        for neighbor in find_neighbors(&current, &pathable, min_x, max_x, min_y, max_y) {
            let tentative_g_score = *g_score.get(&current).unwrap() + 1f32;
            if tentative_g_score < *g_score.get(&neighbor).unwrap_or(&f32::INFINITY) {
                came_from.insert(neighbor, current);
                g_score.insert(neighbor, tentative_g_score);
                let neighbor_f_score = tentative_g_score + h(neighbor, goal);
                open_set.push((Reverse(OrderedFloat::from(neighbor_f_score)), neighbor));
            }
        }
    }
    return None;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unreachable() {
        let start = (0, 0);
        let goal = (2, 2);
        let pathable = vec![
            vec![true, true, false],
            vec![true, false, false],
            vec![false, false, true],
        ];
        let min_x = 0;
        let max_x = 2;
        let min_y = 0;
        let max_y = 2;
        let res = astar(start, goal, &pathable, min_x, max_x, min_y, max_y);
        assert!(res.is_none());
    }

    #[test]
    fn finds_a_path() {
        let start = (0, 0);
        let goal = (2, 2);
        let pathable = vec![
            vec![true, true, false],
            vec![true, false, false],
            vec![true, true, true],
        ];
        let min_x = 0;
        let max_x = 2;
        let min_y = 0;
        let max_y = 2;
        let res = astar(start, goal, &pathable, min_x, max_x, min_y, max_y);
        assert!(res.is_some());
        let res = res.unwrap();
        assert!(res == vec![(0, 0), (0, 1), (1, 2), (2, 2)]);
    }
}
