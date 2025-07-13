use crate::astar::astar;
use crate::navigability_mask::Coords;
use crate::navigability_mask::NavigabilityMask;
use ordered_float::OrderedFloat;
use std::{
    cmp::Reverse,
    collections::{BinaryHeap, HashMap},
};

type TransitionId = u64;

fn join_transition_id(tile: Coords) -> TransitionId {
    ((tile.0 as u64) << 32) | (tile.1 as u64)
}

fn split_transition_id(uid: TransitionId) -> Coords {
    ((uid >> 32) as u32, uid as u32)
}

struct Edge {
    cost: f32,
    destination: TransitionId,
    path: Vec<Coords>,
}

struct Transition {
    position: Coords,
    edges: Vec<Edge>,
}

struct Cluster {
    transitions: Vec<TransitionId>,
}

pub struct HPAStar {
    navigability_mask: NavigabilityMask,
    clusters: Vec<Cluster>, // row-major order
    transitions: HashMap<TransitionId, Transition>,
    width: u32,
    cluster_size: u32,
    entrance_width_breakpoint: u32,
}

impl HPAStar {
    pub fn new(
        navigability_mask: NavigabilityMask,
        cluster_size: u32,
        entrance_width_breakpoint: u32,
    ) -> Self {
        let height = navigability_mask.height as u32;
        assert!(
            height % cluster_size == 0,
            "Grid height was not a multiple of cluster size, was: {}",
            height
        );
        let width = navigability_mask.width as u32;
        assert!(
            width % cluster_size == 0,
            "Grid width was not a multiple of cluster size, was: {}",
            width
        );

        let mut clusters = Vec::<Cluster>::new();
        let mut transitions = HashMap::new();

        // initialize clusters
        for _y in 0..(height / cluster_size) {
            for _x in 0..(width / cluster_size) {
                clusters.push(Cluster {
                    transitions: Vec::new(),
                });
            }
        }

        let insert_transition_pair =
            |a_pos,
             b_pos,
             transitions: &mut HashMap<TransitionId, Transition>,
             clusters: &mut Vec<Cluster>| {
                let insert_empty_transition =
                    |position,
                     transitions: &mut HashMap<TransitionId, Transition>,
                     clusters: &mut Vec<Cluster>| {
                        let cluster_idx = HPAStar::locate_cluster(position, width, cluster_size);
                        let id = join_transition_id((position.0 as u32, position.1 as u32));
                        if transitions.get(&id).is_none() {
                            let edges = vec![];
                            let t = Transition { position, edges };
                            transitions.insert(id, t);
                            clusters[cluster_idx as usize].transitions.push(id);
                        }

                        id
                    };

                let a = insert_empty_transition(a_pos, transitions, clusters);
                let b = insert_empty_transition(b_pos, transitions, clusters);
                transitions.get_mut(&a).unwrap().edges.push(Edge {
                    cost: 1f32,
                    destination: b,
                    path: vec![a_pos, b_pos],
                });
                transitions.get_mut(&b).unwrap().edges.push(Edge {
                    cost: 1f32,
                    destination: a,
                    path: vec![b_pos, a_pos],
                });
            };

        // Iterate over vertical border edges
        for border_start_y in (0..(height - 1)).step_by(cluster_size as usize) {
            // x coord of the border on the right side cluster
            for border_right_x in (cluster_size..(width - 1)).step_by(cluster_size as usize) {
                let border_left_x = border_right_x - 1;
                let mut entrance_start = border_start_y;

                // walk down the edge
                // include 1 extra tile that we'll just treat as if it was closed
                for y in border_start_y..=(border_start_y + cluster_size) {
                    if y < border_start_y + cluster_size
                        && navigability_mask.is_navigable((border_left_x, y))
                        && navigability_mask.is_navigable((border_right_x, y))
                    {
                        // if tile pair is open, keep accumulating
                        continue;
                    }

                    // otherwise close the entrance
                    let first_y = entrance_start;
                    let last_y = y.saturating_sub(1);
                    let entrance_width = y.saturating_sub(first_y);
                    if entrance_width >= entrance_width_breakpoint {
                        let a_pos = (border_left_x, first_y);
                        let b_pos = (border_right_x, first_y);
                        insert_transition_pair(a_pos, b_pos, &mut transitions, &mut clusters);
                        let c_pos = (border_left_x, last_y);
                        let d_pos = (border_right_x, last_y);
                        insert_transition_pair(c_pos, d_pos, &mut transitions, &mut clusters);
                    } else if entrance_width > 0 {
                        let center_y = first_y + ((last_y - first_y) as f32 / 2f32) as u32;
                        let a_pos = (border_left_x, center_y);
                        let b_pos = (border_right_x, center_y);
                        insert_transition_pair(a_pos, b_pos, &mut transitions, &mut clusters);
                    }
                    entrance_start = y + 1;
                }
            }
        }

        // Do the same for all horizontal border edges
        for border_start_x in (0..(width - 1)).step_by(cluster_size as usize) {
            for border_bottom_y in (cluster_size..(height - 1)).step_by(cluster_size as usize) {
                let border_top_y = border_bottom_y - 1;
                let mut entrance_start = border_start_x;

                for x in border_start_x..=(border_start_x + cluster_size) {
                    if x < border_start_x + cluster_size
                        && navigability_mask.is_navigable((x, border_top_y))
                        && navigability_mask.is_navigable((x, border_bottom_y))
                    {
                        continue;
                    }

                    let first_x = entrance_start;
                    let last_x = x.saturating_sub(1);
                    let entrance_width = x.saturating_sub(first_x);
                    if entrance_width >= entrance_width_breakpoint {
                        let a_pos = (first_x, border_top_y);
                        let b_pos = (first_x, border_bottom_y);
                        insert_transition_pair(a_pos, b_pos, &mut transitions, &mut clusters);
                        let c_pos = (last_x, border_top_y);
                        let d_pos = (last_x, border_bottom_y);
                        insert_transition_pair(c_pos, d_pos, &mut transitions, &mut clusters);
                    } else if entrance_width > 0 {
                        let center_x = first_x + ((last_x - first_x) as f32 / 2f32) as u32;
                        let a_pos = (center_x, border_top_y);
                        let b_pos = (center_x, border_bottom_y);
                        insert_transition_pair(a_pos, b_pos, &mut transitions, &mut clusters);
                    }
                    entrance_start = x + 1;
                }
            }
        }

        // for each cluster add intra-edges for each pair of transition nodes
        for cluster_idx in 0..clusters.len() {
            let cluster = &clusters[cluster_idx];
            for i in 0..cluster.transitions.len() {
                for j in (i + 1)..cluster.transitions.len() {
                    let [t1, t2] = {
                        let [maybe_t1, maybe_t2] = transitions
                            .get_disjoint_mut([&cluster.transitions[i], &cluster.transitions[j]]);
                        [maybe_t1.unwrap(), maybe_t2.unwrap()]
                    };
                    let start = t1.position;
                    let goal = t2.position;
                    let pos =
                        HPAStar::get_cluster_position(cluster_idx as u32, width, cluster_size);
                    let min_x = pos.0;
                    let max_x = pos.0 + cluster_size - 1;
                    let min_y = pos.1;
                    let max_y = pos.1 + cluster_size - 1;
                    let maybe_path =
                        astar(start, goal, &navigability_mask, min_x, max_x, min_y, max_y);
                    if let Some(path) = maybe_path {
                        let mut path_rev = path.clone();
                        let cost = (path.len() - 1) as f32;
                        path_rev.reverse();
                        t1.edges.push(Edge {
                            cost,
                            destination: cluster.transitions[j],
                            path,
                        });
                        t2.edges.push(Edge {
                            cost,
                            destination: cluster.transitions[i],
                            path: path_rev,
                        });
                    }
                }
            }
        }

        Self {
            navigability_mask,
            clusters,
            transitions,
            width,
            cluster_size,
            entrance_width_breakpoint,
        }
    }

    fn dist(a: (f32, f32), b: (f32, f32)) -> f32 {
        ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt()
    }

    // Returns the sequence of edges used to get to current in reverse order
    // [nth edge, n-1th edge, .., 1st edge]
    fn reconstruct_path<'a>(
        &self,
        came_from: &HashMap<TransitionId, (TransitionId, &'a Edge)>,
        mut current: TransitionId,
    ) -> Vec<&'a Edge> {
        let mut total_path = vec![];
        while let Some(&prev) = came_from.get(&current) {
            current = prev.0;
            total_path.push(prev.1);
        }
        total_path
    }

    fn h(&self, a_id: TransitionId, b_id: TransitionId) -> f32 {
        let a = self.transitions.get(&a_id).unwrap().position;
        let b = self.transitions.get(&b_id).unwrap().position;
        // (a.0.abs_diff(b.0) + a.1.abs_diff(b.1)) as f32
        ((b.0 as f32 - a.0 as f32).powi(2) + (b.1 as f32 - a.1 as f32).powi(2)).sqrt()
    }

    pub fn find_path(&mut self, start_pos: (u32, u32), end_pos: (u32, u32)) -> Vec<Coords> {
        // check if transition already exists at start/end
        let start_id = self.add_temp_transition(start_pos);
        let end_id = self.add_temp_transition(end_pos);

        // find the path of transitions
        let mut came_from = HashMap::<TransitionId, (TransitionId, &Edge)>::new();

        let mut g_score = HashMap::<TransitionId, f32>::new();
        g_score.insert(start_id, 0f32);

        let mut open_set = BinaryHeap::<(Reverse<OrderedFloat<f32>>, TransitionId)>::new();
        open_set.push((
            Reverse(OrderedFloat::from(self.h(start_id, end_id))),
            start_id,
        ));

        let mut maybe_abstract_path: Option<Vec<&Edge>> = None;
        while !open_set.is_empty() {
            let current = open_set.pop().unwrap().1;
            if current == end_id {
                maybe_abstract_path = Some(self.reconstruct_path(&came_from, current));
                break;
            }

            for neighbor_edge in &self.transitions.get(&current).unwrap().edges {
                let tentative_g_score = *g_score.get(&current).unwrap() + neighbor_edge.cost;
                if tentative_g_score
                    < *g_score
                        .get(&neighbor_edge.destination)
                        .unwrap_or(&f32::INFINITY)
                {
                    came_from.insert(neighbor_edge.destination, (current, neighbor_edge));
                    g_score.insert(neighbor_edge.destination, tentative_g_score);
                    let neighbor_f_score =
                        tentative_g_score + self.h(neighbor_edge.destination, end_id);
                    open_set.push((
                        Reverse(OrderedFloat::from(neighbor_f_score)),
                        neighbor_edge.destination,
                    ));
                }
            }
        }

        // join together all subpaths
        let mut concrete_path: Vec<Coords> = vec![];
        if let Some(mut abstract_path) = maybe_abstract_path {
            while let Some(edge) = abstract_path.pop() {
                let mut clone = edge.path.clone();
                clone.pop(); // avoid double tiles
                concrete_path.append(&mut clone);
            }
        }
        concrete_path.push(end_pos);

        // remove transitions
        self.remove_temp_transition(start_id);
        self.remove_temp_transition(end_id);

        concrete_path
    }

    pub fn add_obstructions(coords: Vec<Coords>) {
        // for each coord, locate_cluster
        // for each unique affected cluster
        //   recalculate transitions
    }

    pub fn remove_obstructions(coords: Vec<Coords>) {
        // for each coord, locate_cluster
        // for each unique affected cluster
        //   recalculate transitions
    }

    fn add_temp_transition(&mut self, p: (u32, u32)) -> TransitionId {
        let id = join_transition_id(p);
        if self.transitions.get(&id).is_some() {
            return id;
        }
        let edges = vec![];
        let transition = Transition { edges, position: p };
        self.transitions.insert(id, transition);

        let cluster_idx = HPAStar::locate_cluster(p, self.width, self.cluster_size) as usize;
        let cluster = &self.clusters[cluster_idx];
        for i in 0..cluster.transitions.len() {
            let [t1, t2] = {
                let [maybe_t1, maybe_t2] = self
                    .transitions
                    .get_disjoint_mut([&cluster.transitions[i], &id]);
                [maybe_t1.unwrap(), maybe_t2.unwrap()]
            };
            let start = t1.position;
            let goal = t2.position;
            let pos =
                HPAStar::get_cluster_position(cluster_idx as u32, self.width, self.cluster_size);
            let min_x = pos.0;
            let max_x = pos.0 + self.cluster_size - 1;
            let min_y = pos.1;
            let max_y = pos.1 + self.cluster_size - 1;
            let maybe_path = astar(
                start,
                goal,
                &self.navigability_mask,
                min_x,
                max_x,
                min_y,
                max_y,
            );
            if let Some(path) = maybe_path {
                let mut path_rev = path.clone();
                let cost = (path.len() - 1) as f32;
                path_rev.reverse();
                t1.edges.push(Edge {
                    cost,
                    destination: id,
                    path,
                });
                t2.edges.push(Edge {
                    cost,
                    destination: cluster.transitions[i],
                    path: path_rev,
                });
            }
        }
        return id;
    }

    fn remove_temp_transition(&mut self, removee_id: TransitionId) {
        let removee = self.transitions.get(&removee_id).unwrap();
        let cluster_idx = HPAStar::locate_cluster(removee.position, self.width, self.cluster_size);
        let cluster_transitions = &self.clusters[cluster_idx as usize].transitions;
        for t_id in cluster_transitions {
            let t = self.transitions.get_mut(t_id).unwrap();
            if let Some(pos) = t
                .edges
                .iter()
                .position(|&Edge { destination, .. }| destination == removee_id)
            {
                t.edges.remove(pos);
            }
        }
        self.transitions.remove(&removee_id);
    }

    fn locate_cluster(p: (u32, u32), width: u32, cluster_size: u32) -> u32 {
        let x = (p.0 as f32 / cluster_size as f32).floor() as u32;
        let y = (p.1 as f32 / cluster_size as f32).floor() as u32;
        y * (width / cluster_size) + x
    }

    fn get_cluster_position(idx: u32, width: u32, cluster_size: u32) -> Coords {
        let clusters_per_row = width / cluster_size;
        let x = idx % clusters_per_row;
        let y = idx / clusters_per_row;
        (x * cluster_size, y * cluster_size)
    }
}

#[cfg(test)]
mod test_helpers {
    use crate::navigability_mask::Coords;
    use plotters::prelude::*;
    use std::collections::HashSet;
    use std::fs;
    use std::fs::File;
    use std::io::{self, BufRead};
    use std::path::Path;

    pub fn load_map<P: AsRef<Path>>(path: P) -> io::Result<Vec<Vec<bool>>> {
        let file = File::open(path)?;
        let mut lines = io::BufReader::new(file).lines();

        // Read headers
        let mut map_type = String::new();
        let mut width = 0;
        let mut height = 0;
        while let Some(Ok(line)) = lines.next() {
            if line.starts_with("type ") {
                map_type = line["type ".len()..].to_string();
            } else if line.starts_with("height ") {
                height = line["height ".len()..].trim().parse().unwrap();
            } else if line.starts_with("width ") {
                width = line["width ".len()..].trim().parse().unwrap();
            } else if line.trim() == "map" {
                break;
            }
        }
        assert_eq!(map_type, "octile", "Unsupported map type");

        // Read rows
        let mut grid = Vec::with_capacity(height);
        for _ in 0..height {
            let row = lines.next().unwrap().unwrap();
            assert_eq!(row.len(), width);
            let passable_row = row
                .chars()
                .map(|c| match c {
                    '.' | 'G' | 'S' => true,        // passable, swamp allowed
                    'T' | 'W' | '@' | 'O' => false, // trees, water, out-of-bounds
                    _ => false,
                })
                .collect();
            grid.push(passable_row);
        }

        Ok(grid)
    }

    pub fn draw_grid(
        grid: Vec<Vec<bool>>,
        path: Vec<Coords>,
        transitions: Vec<Coords>,
        cluster_size: u32,
    ) {
        let width = grid[0].len() as u32;
        let height = grid.len() as u32;
        let pixel_size = 4;
        let img_width = width * pixel_size;
        let img_height = height * pixel_size;

        fs::create_dir_all("test-output").unwrap();
        let backend = BitMapBackend::new("test-output/grid.png", (img_width, img_height));
        let root = backend.into_drawing_area();
        root.fill(&WHITE).unwrap();

        let transition_set: HashSet<_> = transitions.into_iter().collect();
        let path_set: HashSet<_> = path.into_iter().collect();

        for (y, row) in grid.iter().enumerate() {
            for (x, &walkable) in row.iter().enumerate() {
                let x = x as u32;
                let y = y as u32;
                let x0 = x * pixel_size;
                let y0 = y * pixel_size;
                let x1 = x0 + pixel_size;
                let y1 = y0 + pixel_size;

                let mut color = WHITE;
                if !walkable {
                    color = BLACK;
                } else if transition_set.contains(&(x, y)) {
                    color = RGBColor(150, 180, 210); // light blue
                } else if path_set.contains(&(x, y)) {
                    color = RED;
                }

                root.draw(&Rectangle::new(
                    [(x0 as i32, y0 as i32), (x1 as i32, y1 as i32)],
                    color.filled(),
                ))
                .unwrap();
            }
        }

        // Draw vertical grid lines
        for x in 0..=width {
            let xpos = (x * pixel_size) as i32;
            root.draw(&PathElement::new(
                vec![(xpos, 0), (xpos, img_height as i32)],
                &RGBAColor(0, 0, 0, 0.2),
            ))
            .unwrap();
        }

        // Draw horizontal grid lines
        for y in 0..=height {
            let ypos = (y * pixel_size) as i32;
            root.draw(&PathElement::new(
                vec![(0, ypos), (img_width as i32, ypos)],
                &RGBAColor(0, 0, 0, 0.2),
            ))
            .unwrap();
        }

        // Cluster grid lines
        let thick_style = ShapeStyle::from(&BLUE).stroke_width(1);
        for x in (0..=width).step_by(cluster_size as usize) {
            let xpos = (x * pixel_size) as i32;
            root.draw(&PathElement::new(
                vec![(xpos, 0), (xpos, img_height as i32)],
                thick_style.clone(),
            ))
            .unwrap();
        }
        for y in (0..=height).step_by(cluster_size as usize) {
            let ypos = (y * pixel_size) as i32;
            root.draw(&PathElement::new(
                vec![(0, ypos), (img_width as i32, ypos)],
                thick_style.clone(),
            ))
            .unwrap();
        }
    }

    pub fn generate_test_grid() -> Vec<Vec<bool>> {
        let width = 512;
        let height = 512;

        let mut grid = vec![vec![true; width]; height];

        // Add some barriers: vertical walls every 64 columns, leaving gaps
        for x in (64..width).step_by(64) {
            for y in 0..height {
                if y % 20 != 0 {
                    // leave periodic gaps
                    grid[y][x] = false;
                }
            }
        }

        // Add horizontal walls every 64 rows, leaving gaps
        for y in (64..height).step_by(64) {
            for x in 0..width {
                if x % 20 != 0 {
                    grid[y][x] = false;
                }
            }
        }

        grid
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn one_cluster() {
        let mut pathable = vec![];
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        for _j in 0..cluster_size {
            let mut row = vec![];
            for _i in 0..cluster_size {
                row.push(true);
            }
            pathable.push(row);
        }
        let nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable),
            cluster_size,
            entrance_width_breakpoint,
        );
        assert!(nav.clusters.len() == 1);
        assert!(nav.transitions.len() == 0);
    }

    #[test]
    fn two_clusters() {
        let mut pathable = vec![];
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        for _j in 0..cluster_size {
            let mut row = vec![];
            for _i in 0..cluster_size * 2 {
                row.push(true);
            }
            pathable.push(row);
        }
        let nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable),
            cluster_size,
            entrance_width_breakpoint,
        );
        assert!(nav.clusters.len() == 2);
        assert!(
            nav.transitions.len() == 4,
            "Expected 4 transitions, was: {}",
            nav.transitions.len()
        );
        for t in nav.transitions.iter() {
            assert!(t.1.edges.len() == 2);
        }
    }

    #[test]
    fn four_clusters() {
        let mut pathable = vec![];
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        for _j in 0..cluster_size * 2 {
            let mut row = vec![];
            for _i in 0..cluster_size * 2 {
                row.push(true);
            }
            pathable.push(row);
        }
        let nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable),
            cluster_size,
            entrance_width_breakpoint,
        );
        assert!(nav.clusters.len() == 4);
        assert!(
            nav.transitions.len() == 12,
            "Expected 12 transitions, was: {}",
            nav.transitions.len()
        );
    }

    #[test]
    fn four_clusters_with_obstacles() {
        let mut pathable = vec![];
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        for _j in 0..cluster_size * 2 {
            let mut row = vec![];
            for _i in 0..cluster_size * 2 {
                row.push(true);
            }
            pathable.push(row);
        }

        // left wall unpathable
        for y in 0..cluster_size * 2 {
            pathable[y as usize][0] = false;
        }

        // top vertical wall closed
        for y in 0..cluster_size {
            pathable[y as usize][cluster_size as usize] = false;
        }

        // left horizontal wall semi-closed
        for x in 4..cluster_size {
            pathable[cluster_size as usize][x as usize] = false;
        }

        // right horizontal wall split
        for x in cluster_size + 4..cluster_size * 2 - 4 {
            pathable[cluster_size as usize][x as usize] = false;
        }

        let mut nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable.clone()),
            cluster_size,
            entrance_width_breakpoint,
        );
        assert!(
            nav.clusters[0].transitions.len() == 1,
            "Expected 1 transition, was: {}",
            nav.clusters[0].transitions.len(),
        );
        assert!(
            nav.clusters[1].transitions.len() == 2,
            "Expected 2 transitions, was: {}",
            nav.clusters[1].transitions.len(),
        );
        assert!(
            nav.clusters[2].transitions.len() == 3,
            "Expected 4 transitions, was: {}",
            nav.clusters[2].transitions.len(),
        );
        assert!(
            nav.clusters[3].transitions.len() == 4,
            "Expected 4 transitions, was: {}",
            nav.clusters[3].transitions.len(),
        );

        assert!(
            nav.clusters[0]
                .transitions
                .iter()
                .all(|t_id| { nav.transitions.get(t_id).unwrap().edges.len() == 1 }),
            "Expected all transitions in cluster 0 to have 1 edge"
        );
        assert!(
            nav.clusters[1]
                .transitions
                .iter()
                .all(|t_id| { nav.transitions.get(t_id).unwrap().edges.len() == 2 }),
            "Expected all transitions in cluster 1 to have 2 edges"
        );
        assert!(
            nav.clusters[2]
                .transitions
                .iter()
                .all(|t_id| { nav.transitions.get(t_id).unwrap().edges.len() == 3 }),
            "Expected all transitions in cluster 2 to have 3 edges"
        );
        assert!(
            nav.clusters[3]
                .transitions
                .iter()
                .all(|t_id| { nav.transitions.get(t_id).unwrap().edges.len() == 4 }),
            "Expected all transitions in cluster 3 to have 4 edges"
        );

        let start: Coords = (5, 5);
        let end: Coords = (cluster_size + 5, 5);
        let hpastar_path = nav.find_path(start, end);
        /* let astar_path = astar::astar(
            start,
            end,
            &pathable,
            0,
            CLUSTER_SIZE * 2 - 1,
            0,
            CLUSTER_SIZE * 2 - 1,
        )
        .unwrap();
        */

        let mut seen = HashSet::<Coords>::new();
        for pos in &hpastar_path {
            if !seen.insert(*pos) {
                panic!("The path visits {:?} more than once", pos);
            }
        }

        test_helpers::draw_grid(
            pathable,
            hpastar_path,
            nav.transitions
                .into_iter()
                .map(|(_id, transition)| transition.position)
                .collect(),
            cluster_size,
        );
    }

    #[test]
    #[ignore]
    fn big_grid() {
        let pathable = test_helpers::generate_test_grid();
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        let mut nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable.clone()),
            cluster_size,
            entrance_width_breakpoint,
        );
        let start: (u32, u32) = (8, 8);
        let end: (u32, u32) = (500, 500);
        let hpastar_path = nav.find_path(start, end);

        assert!(
            hpastar_path.len() > 500,
            "path length was: {}",
            hpastar_path.len()
        );
        test_helpers::draw_grid(
            pathable,
            hpastar_path,
            nav.transitions
                .into_iter()
                .map(|(_id, transition)| transition.position)
                .collect(),
            cluster_size,
        );
    }

    #[test]
    #[ignore]
    fn real_map() -> Result<(), Box<dyn std::error::Error>> {
        let path = "benchmarks/BigGameHunters.map";
        let pathable = test_helpers::load_map(path)?;
        let cluster_size = 16;
        let entrance_width_breakpoint = 8;
        let mut nav = HPAStar::new(
            NavigabilityMask::from_row_major_vec(pathable.clone()),
            cluster_size,
            entrance_width_breakpoint,
        );
        let start: (u32, u32) = (20, 16);
        let end: (u32, u32) = (132, 507);
        let hpastar_path = nav.find_path(start, end);
        test_helpers::draw_grid(
            pathable,
            hpastar_path,
            nav.transitions
                .into_iter()
                .map(|(_id, transition)| transition.position)
                .collect(),
            cluster_size,
        );
        Ok(())
    }
}
