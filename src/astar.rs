use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
    u32,
};

use crate::navigability_mask::{self, Coords, NavigabilityMask};

// includes both start and destination elements
fn reconstruct_path(came_from: &HashMap<Coords, Coords>, mut current: Coords) -> Vec<Coords> {
    let mut total_path = vec![current];
    while let Some(&prev) = came_from.get(&current) {
        current = prev;
        total_path.push(current);
    }
    total_path.reverse();
    total_path
}

// Octile distance
const SCALE: u32 = 1000;
pub const STRAIGHT_COST: u32 = 1000;
pub const DIAGONAL_COST: u32 = 1414; // approx sqrt(2) * 1000

pub fn h(a: Coords, b: Coords) -> u32 {
    let dx = (b.0 as i32 - a.0 as i32).abs() as u32;
    let dy = (b.1 as i32 - a.1 as i32).abs() as u32;
    let min = dx.min(dy);
    let max = dx.max(dy);
    DIAGONAL_COST * min + STRAIGHT_COST * (max - min)
}

pub fn astar(
    start: Coords,
    goal: Coords,
    navigability_mask: &NavigabilityMask,
    min_x: u32,
    max_x: u32,
    min_y: u32,
    max_y: u32,
) -> Option<Vec<Coords>> {
    let mut came_from = HashMap::<Coords, Coords>::new();

    let mut g_score = HashMap::<Coords, u32>::new();
    g_score.insert(start, 0u32);

    let mut open_set = BinaryHeap::<(Reverse<u32>, Coords)>::new();
    open_set.push((Reverse(h(start, goal)), start));

    while !open_set.is_empty() {
        let current = open_set.pop().unwrap().1;
        if current == goal {
            return Some(reconstruct_path(&came_from, current));
        }

        for (neighbor, is_diag) in navigability_mask.get_navigable_neighbors(
            current,
            min_x as usize,
            max_x as usize,
            min_y as usize,
            max_y as usize,
        ) {
            let tentative_g_score = *g_score.get(&current).unwrap()
                + if is_diag {
                    DIAGONAL_COST
                } else {
                    STRAIGHT_COST
                };
            if tentative_g_score < *g_score.get(&neighbor).unwrap_or(&u32::MAX) {
                came_from.insert(neighbor, current);
                g_score.insert(neighbor, tentative_g_score);
                let neighbor_f_score = tentative_g_score + h(neighbor, goal);
                open_set.push((Reverse(neighbor_f_score), neighbor));
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
        let pathable = NavigabilityMask::from_row_major_vec(vec![
            vec![true, true, false],
            vec![true, false, false],
            vec![false, false, true],
        ]);
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
        let pathable = NavigabilityMask::from_row_major_vec(vec![
            vec![true, true, false],
            vec![true, false, false],
            vec![true, true, true],
        ]);
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
