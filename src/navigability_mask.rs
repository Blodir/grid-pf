use bitvector::BitVector;

pub type Coords = (u32, u32);

pub struct NavigabilityMask {
    bitvec: BitVector,
    pub height: usize,
    pub width: usize,
}
impl NavigabilityMask {
    pub fn from_row_major_vec(pathable: Vec<Vec<bool>>) -> Self {
        let height = pathable.len() as usize;
        let width = pathable[0].len() as usize;
        let mut bitvec = BitVector::new(height * width);
        for (row_idx, row) in pathable.iter().enumerate() {
            for (col_idx, elem) in row.iter().enumerate() {
                if *elem {
                    bitvec.insert(row_idx * width + col_idx);
                }
            }
        }
        Self {
            bitvec,
            height,
            width,
        }
    }
    pub fn is_navigable(&self, (x, y): Coords) -> bool {
        self.bitvec.contains(y as usize * self.width + x as usize)
    }

    pub fn get_navigable_neighbors<'a>(
        &'a self,
        pos: Coords,
        min_x: usize,
        max_x: usize,
        min_y: usize,
        max_y: usize,
    ) -> impl Iterator<Item = Coords> + 'a {
        let (x, y) = (pos.0 as i32, pos.1 as i32);

        const DELTAS: [(i32, i32); 8] = [
            (-1, 0),
            (1, 0),
            (0, -1),
            (0, 1),
            (-1, -1),
            (1, -1),
            (-1, 1),
            (1, 1),
        ];

        DELTAS.into_iter().filter_map(move |(dx, dy)| {
            let nx = x + dx;
            let ny = y + dy;

            if nx >= min_x as i32 && nx <= max_x as i32 && ny >= min_y as i32 && ny <= max_y as i32
            {
                let p = (nx as u32, ny as u32);
                if self.is_navigable(p) {
                    return Some(p);
                }
            }
            None
        })
    }

    pub fn set_rect(&mut self, top_left: Coords, bot_right: Coords) {}
}
