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

    /// returns (neighbor position, neighbor is diagonal)
    pub fn get_navigable_neighbors<'a>(
        &'a self,
        pos: Coords,
        min_x: usize,
        max_x: usize,
        min_y: usize,
        max_y: usize,
    ) -> impl Iterator<Item = (Coords, bool)> + 'a {
        let (x, y) = (pos.0 as i32, pos.1 as i32);

        const DELTAS: [((i32, i32), bool); 8] = [
            ((-1, 0), false),
            ((1, 0), false),
            ((0, -1), false),
            ((0, 1), false),
            ((-1, -1), true),
            ((1, -1), true),
            ((-1, 1), true),
            ((1, 1), true),
        ];

        DELTAS.into_iter().filter_map(move |((dx, dy), is_diag)| {
            let nx = x + dx;
            let ny = y + dy;

            if nx >= min_x as i32 && nx <= max_x as i32 && ny >= min_y as i32 && ny <= max_y as i32
            {
                let p = (nx as u32, ny as u32);
                if self.is_navigable(p) {
                    return Some((p, is_diag));
                }
            }

            None
        })
    }

    pub fn set_rect(&mut self, top_left: Coords, bot_right: Coords) {}

    // Using Bresenham line drawing algorithm
    pub fn visibility_check(&self, a: Coords, b: Coords) -> bool {
        let x0 = a.0 as i32;
        let y0 = a.1 as i32;
        let x1 = b.0 as i32;
        let y1 = b.1 as i32;

        let dx = (x1 - x0).abs();
        let dy = -(y1 - y0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;

        let mut x = x0;
        let mut y = y0;

        loop {
            if x < 0 || y < 0 || !self.is_navigable((x as u32, y as u32)) {
                return false;
            }
            if x == x1 && y == y1 {
                return true;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x += sx;
            }
            if e2 <= dx {
                err += dx;
                y += sy;
            }
        }
    }
}
