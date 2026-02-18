use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::{Circle, Edge};

/// Common trait for all geometric primitives that can be directly queried in the quadtree
/// for collisions with the edges of the registered hazards. These include: [Rect], [Edge] and [Circle].
pub trait QTQueryable: CollidesWith<Edge> + CollidesWith<Rect> {
    /// Checks which quadrants the entity collides with.
    fn collides_with_quadrants(&self, _r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        debug_assert!(_r.quadrants().iter().zip(qs.iter()).all(|(q, r)| *q == **r));
        qs.map(|q| self.collides_with(q))
    }
}

impl QTQueryable for Circle {}
impl QTQueryable for Rect {}

impl QTQueryable for Edge {
    fn collides_with_quadrants(&self, r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        debug_assert!(r.quadrants().iter().zip(qs.iter()).all(|(q, r)| *q == **r));

        let x0 = self.start.0;
        let y0 = self.start.1;
        let dx = self.end.0 - x0;
        let dy = self.end.1 - y0;
        let dxy = dx * dy;

        if dxy.abs() < 1e-12 {
            // Axis-aligned (or degenerate) edge — equivalent to its bbox
            let (x_min, x_max) = if dx >= 0.0 { (x0, x0 + dx) } else { (x0 + dx, x0) };
            let (y_min, y_max) = if dy >= 0.0 { (y0, y0 + dy) } else { (y0 + dy, y0) };

            if x_min.max(r.x_min) > x_max.min(r.x_max)
                || y_min.max(r.y_min) > y_max.min(r.y_max)
            {
                return [false; 4];
            }

            let cx = (r.x_min + r.x_max) * 0.5;
            let cy = (r.y_min + r.y_max) * 0.5;
            let left = x_min <= cx;
            let right = x_max >= cx;
            let below = y_min <= cy;
            let above = y_max >= cy;
            return [right && above, left && above, left && below, right && below];
        }

        let dx_min_dy = (r.x_min - x0) * dy;
        let dx_max_dy = (r.x_max - x0) * dy;
        let dy_min_dx = (r.y_min - y0) * dx;
        let dy_max_dx = (r.y_max - y0) * dx;

        let r_min = dxy.min(0.0).max(dx_min_dy.min(dx_max_dy)).max(dy_min_dx.min(dy_max_dx));
        let r_max = dxy.max(0.0).min(dx_min_dy.max(dx_max_dy)).min(dy_min_dx.max(dy_max_dx));

        if r_max < r_min {
            return [false; 4];
        }

        let dcx_dy = dx_min_dy + dx_max_dy;
        let dcy_dx = dy_min_dx + dy_max_dx;
        let r_min2 = r_min + r_min;
        let r_max2 = r_max + r_max;

        let q = [
            dcx_dy.max(dcy_dx) <= r_max2,
            r_min2.max(dcy_dx) <= r_max2.min(dcx_dy),
            r_min2 <= dcx_dy.min(dcy_dx),
            r_min2.max(dcx_dy) <= r_max2.min(dcy_dx),
        ];

        match (dx > 0.0, dy > 0.0) {
            (true, true) => [q[0], q[1], q[2], q[3]],
            (false, true) => [q[3], q[2], q[1], q[0]],
            (true, false) => [q[1], q[0], q[3], q[2]],
            (false, false) => [q[2], q[3], q[0], q[1]],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::primitives::{Edge, Point, Rect};
    use std::cmp::Ordering;

    // ==================== OLD IMPLEMENTATION ====================

    const OLD_Q0: u8 = 0b0001;
    const OLD_Q1: u8 = 0b0010;
    const OLD_Q2: u8 = 0b0100;
    const OLD_Q3: u8 = 0b1000;
    const OLD_ALL: u8 = 0b1111;

    fn old_collides_with_quadrants(edge: &Edge, r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        let mut determined: u8 = 0;
        let mut collides: u8 = 0;

        let e_x_min = edge.x_min();
        let e_x_max = edge.x_max();
        let e_y_min = edge.y_min();
        let e_y_max = edge.y_max();

        for (i, q) in qs.iter().enumerate() {
            let bit = 1 << i;
            let x_no_overlap = e_x_min.max(q.x_min) > e_x_max.min(q.x_max);
            let y_no_overlap = e_y_min.max(q.y_min) > e_y_max.min(q.y_max);

            if x_no_overlap || y_no_overlap {
                determined |= bit;
            } else if q.collides_with(&edge.start) || q.collides_with(&edge.end) {
                determined |= bit;
                collides |= bit;
            }
        }

        if determined == OLD_ALL {
            return old_bits_to_array(collides);
        }

        let c = r.centroid();
        let [top, left, bottom, right] = r.sides();

        let h_bisect = Edge {
            start: Point(r.x_min, c.1),
            end: Point(r.x_max, c.1),
        };
        let v_bisect = Edge {
            start: Point(c.0, r.y_min),
            end: Point(c.0, r.y_max),
        };

        old_half_intersect_bits(edge, &left, OLD_Q1, OLD_Q2, &mut determined, &mut collides);
        old_half_intersect_bits(edge, &right, OLD_Q3, OLD_Q0, &mut determined, &mut collides);
        old_half_intersect_bits(edge, &top, OLD_Q0, OLD_Q1, &mut determined, &mut collides);
        old_half_intersect_bits(edge, &bottom, OLD_Q2, OLD_Q3, &mut determined, &mut collides);
        old_half_intersect_bits(
            edge,
            &h_bisect,
            OLD_Q1 | OLD_Q2,
            OLD_Q0 | OLD_Q3,
            &mut determined,
            &mut collides,
        );
        old_half_intersect_bits(
            edge,
            &v_bisect,
            OLD_Q2 | OLD_Q3,
            OLD_Q0 | OLD_Q1,
            &mut determined,
            &mut collides,
        );

        old_bits_to_array(collides)
    }

    fn old_bits_to_array(bits: u8) -> [bool; 4] {
        [
            bits & OLD_Q0 != 0,
            bits & OLD_Q1 != 0,
            bits & OLD_Q2 != 0,
            bits & OLD_Q3 != 0,
        ]
    }

    fn old_half_intersect_bits(
        e1: &Edge,
        e2: &Edge,
        fst_mask: u8,
        sec_mask: u8,
        determined: &mut u8,
        collides: &mut u8,
    ) {
        let relevant = (fst_mask | sec_mask) & !*determined;
        if relevant != 0 {
            if let Some((_, e2_col_loc)) = old_edge_intersection_half(e1, e2) {
                let hit_mask = match e2_col_loc {
                    OldCollisionHalf::FirstHalf => fst_mask,
                    OldCollisionHalf::SecondHalf => sec_mask,
                    OldCollisionHalf::Halfway => fst_mask | sec_mask,
                };
                *determined |= hit_mask;
                *collides |= hit_mask;
            }
        }
    }

    fn old_edge_intersection_half(
        e1: &Edge,
        e2: &Edge,
    ) -> Option<(OldCollisionHalf, OldCollisionHalf)> {
        let Point(x1, y1) = e1.start;
        let Point(x2, y2) = e1.end;
        let Point(x3, y3) = e2.start;
        let Point(x4, y4) = e2.end;

        let t_nom = (x2 - x4) * (y4 - y3) - (y2 - y4) * (x4 - x3);
        let t_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
        let u_nom = (x2 - x4) * (y2 - y1) - (y2 - y4) * (x2 - x1);
        let u_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
        if t_denom == 0.0 || u_denom == 0.0 {
            return None;
        }

        let t = t_nom / t_denom;
        let u = u_nom / u_denom;
        if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
            let e1_loc = match t.partial_cmp(&0.5).unwrap() {
                Ordering::Greater => OldCollisionHalf::FirstHalf,
                Ordering::Less => OldCollisionHalf::SecondHalf,
                Ordering::Equal => OldCollisionHalf::Halfway,
            };
            let e2_loc = match u.partial_cmp(&0.5).unwrap() {
                Ordering::Greater => OldCollisionHalf::FirstHalf,
                Ordering::Less => OldCollisionHalf::SecondHalf,
                Ordering::Equal => OldCollisionHalf::Halfway,
            };
            return Some((e1_loc, e2_loc));
        }
        None
    }

    enum OldCollisionHalf {
        FirstHalf,
        Halfway,
        SecondHalf,
    }

    // ==================== TEST HELPERS ====================

    fn make_rect(x_min: f32, y_min: f32, x_max: f32, y_max: f32) -> Rect {
        Rect { x_min, y_min, x_max, y_max }
    }

    fn make_edge(x1: f32, y1: f32, x2: f32, y2: f32) -> Edge {
        Edge { start: Point(x1, y1), end: Point(x2, y2) }
    }

    /// Compare new vs old. Panics on any difference.
    fn compare(edge: &Edge, r: &Rect) {
        let qs = r.quadrants();
        let qs_ref = [&qs[0], &qs[1], &qs[2], &qs[3]];

        let old = old_collides_with_quadrants(edge, r, qs_ref);
        let new = edge.collides_with_quadrants(r, qs_ref);

        if old != new {
            let naive = qs_ref.map(|q| q.collides_with(edge));
            panic!(
                "MISMATCH!\n  edge:  {:?}\n  rect:  {:?}\n  old:   {:?}\n  new:   {:?}\n  naive: {:?}",
                edge, r, old, new, naive
            );
        }
    }

    // ==================== TESTS ====================

    #[test]
    fn test_edge_fully_inside_one_quadrant() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(6.0, 6.0, 9.0, 9.0), &r);
        compare(&make_edge(1.0, 6.0, 4.0, 9.0), &r);
        compare(&make_edge(1.0, 1.0, 4.0, 4.0), &r);
        compare(&make_edge(6.0, 1.0, 9.0, 4.0), &r);
    }

    #[test]
    fn test_edge_crossing_vertical_midline() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(3.0, 7.0, 7.0, 7.0), &r);
        compare(&make_edge(3.0, 3.0, 7.0, 3.0), &r);
    }

    #[test]
    fn test_edge_crossing_horizontal_midline() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(7.0, 3.0, 7.0, 7.0), &r);
        compare(&make_edge(3.0, 3.0, 3.0, 7.0), &r);
    }

    #[test]
    fn test_edge_diagonal_q0_q2() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(1.0, 1.0, 9.0, 9.0), &r);
        compare(&make_edge(1.0, 4.0, 9.0, 6.0), &r);
        compare(&make_edge(2.0, 2.0, 8.0, 8.0), &r);
    }

    #[test]
    fn test_edge_diagonal_q1_q3() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(1.0, 9.0, 9.0, 1.0), &r);
        compare(&make_edge(9.0, 1.0, 1.0, 9.0), &r);
    }

    #[test]
    fn test_edge_on_midlines() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(5.0, 2.0, 5.0, 8.0), &r);
        compare(&make_edge(2.0, 5.0, 8.0, 5.0), &r);
        compare(&make_edge(5.0, 5.0, 8.0, 8.0), &r);
        compare(&make_edge(2.0, 3.0, 5.0, 5.0), &r);
    }

    #[test]
    fn test_edge_outside_rect() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(11.0, 3.0, 15.0, 7.0), &r);
        compare(&make_edge(3.0, 11.0, 7.0, 15.0), &r);
        compare(&make_edge(3.0, -5.0, 7.0, -1.0), &r);
    }

    #[test]
    fn test_edge_partially_outside_rect() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(-5.0, 7.0, 7.0, 7.0), &r);
        compare(&make_edge(3.0, 3.0, 15.0, 3.0), &r);
        compare(&make_edge(-5.0, 5.0, 15.0, 5.0), &r);
        compare(&make_edge(-2.0, -2.0, 12.0, 12.0), &r);
    }

    #[test]
    fn test_edge_touching_rect_border() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(2.0, 10.0, 8.0, 10.0), &r);
        compare(&make_edge(0.0, 2.0, 0.0, 8.0), &r);
        compare(&make_edge(0.0, 0.0, 5.0, 5.0), &r);
    }

    #[test]
    fn test_edge_spanning_all_quadrants() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(1.0, 5.0, 9.0, 5.0), &r);
        compare(&make_edge(5.0, 1.0, 5.0, 9.0), &r);
    }

    #[test]
    fn test_non_square_rect() {
        let r = make_rect(0.0, 0.0, 20.0, 10.0);
        compare(&make_edge(3.0, 3.0, 17.0, 7.0), &r);
        compare(&make_edge(8.0, 3.0, 12.0, 7.0), &r);

        let r = make_rect(0.0, 0.0, 10.0, 20.0);
        compare(&make_edge(3.0, 3.0, 7.0, 17.0), &r);
        compare(&make_edge(3.0, 8.0, 7.0, 12.0), &r);
    }

    #[test]
    fn test_edge_along_quadrant_boundary() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(5.0, 6.0, 5.0, 9.0), &r);
        compare(&make_edge(5.0, 1.0, 5.0, 4.0), &r);
        compare(&make_edge(6.0, 5.0, 9.0, 5.0), &r);
        compare(&make_edge(1.0, 5.0, 4.0, 5.0), &r);
    }

    #[test]
    fn test_tiny_edges() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        compare(&make_edge(4.99, 4.99, 5.01, 5.01), &r);
        compare(&make_edge(4.99, 7.0, 5.01, 7.0), &r);
    }

    #[test]
    fn test_edge_bbox_spans_both_midlines_but_hits_few_quadrants() {
        let r = make_rect(0.0, 0.0, 10.0, 10.0);
        // bbox spans both midlines (cx=5,cy=5) but edge only hits 1 quadrant
        compare(&make_edge(-2.0, 3.0, 5.0, 12.0), &r);   // Q1 only
        compare(&make_edge(-2.0, 5.0, 5.0, -2.0), &r);    // Q2 only
        compare(&make_edge(-2.0, 5.0, 6.0, -2.0), &r);    // Q2 only
        compare(&make_edge(-2.0, 4.0, 7.0, 12.0), &r);    // Q1 only
    }

    #[test]
    fn test_stress_random_edges() {
        let r = make_rect(-10.0, -10.0, 10.0, 10.0);
        let coords: Vec<f32> = vec![
            -15.0, -10.0, -7.5, -5.0, -2.5, 0.0, 2.5, 5.0, 7.5, 10.0, 15.0,
        ];

        let mut count = 0;
        for &x1 in &coords {
            for &y1 in &coords {
                for &x2 in &coords {
                    for &y2 in &coords {
                        if x1 == x2 && y1 == y2 {
                            continue;
                        }
                        compare(&make_edge(x1, y1, x2, y2), &r);
                        count += 1;
                    }
                }
            }
        }
        eprintln!("Tested {count} edge/rect combinations");
    }

    #[test]
    fn test_stress_non_square_random() {
        let r = make_rect(-5.0, -10.0, 15.0, 10.0);
        let coords: Vec<f32> = vec![-12.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0];

        let mut count = 0;
        for &x1 in &coords {
            for &y1 in &coords {
                for &x2 in &coords {
                    for &y2 in &coords {
                        if x1 == x2 && y1 == y2 {
                            continue;
                        }
                        compare(&make_edge(x1, y1, x2, y2), &r);
                        count += 1;
                    }
                }
            }
        }
        eprintln!("Tested {count} edge/rect combinations (non-square)");
    }
}
