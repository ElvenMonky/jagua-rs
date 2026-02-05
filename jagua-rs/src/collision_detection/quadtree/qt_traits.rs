use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::{Circle, Edge, Point, SPolygon};
use std::cmp::Ordering;
use std::f32::consts::PI;

/// Common trait for all geometric primitives that can be directly queried in the quadtree
/// for collisions with the edges of the registered hazards. These include: [Rect], [Edge], [Circle], and [SPolygon].
pub trait QTQueryable: CollidesWith<Edge> + CollidesWith<Rect> {
    /// Checks which quadrants the entity collides with.
    fn collides_with_quadrants(&self, _r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        debug_assert!(_r.quadrants().iter().zip(qs.iter()).all(|(q, r)| *q == **r));
        qs.map(|q| self.collides_with(q))
    }

    /// Returns true if pigeonhole principle guarantees collision:
    /// entity_presence_area + haz_presence_area > bbox.area()
    fn guarantees_collision(&self, _bbox: &Rect, _haz_presence_area: f32) -> bool {
        false
    }
}

#[inline]
fn circle_rect_intersection_area(cx: f32, cy: f32, r: f32, rect: &Rect) -> f32 {
    let r_sq = r * r;

    // Fast negative: circle doesn't touch rect
    let cx_clamped = cx.clamp(rect.x_min, rect.x_max);
    let cy_clamped = cy.clamp(rect.y_min, rect.y_max);
    let d_sq = (cx - cx_clamped).powi(2) + (cy - cy_clamped).powi(2);
    if d_sq > r_sq {
        return 0.0;
    }

    if cx - r >= rect.x_min && cx + r <= rect.x_max && cy - r >= rect.y_min && cy + r <= rect.y_max {
        return PI * r_sq;
    }

    // Intersection interval along x (clamped by rect)
    let dy = (rect.y_min - cy).max(cy - rect.y_max).max(0.0);
    let half_x = (r_sq - dy * dy).max(0.0).sqrt();
    let ix_min = (cx - half_x).max(rect.x_min);
    let ix_max = (cx + half_x).min(rect.x_max);

    // Intersection interval along y (clamped by rect)
    let dx = (rect.x_min - cx).max(cx - rect.x_max).max(0.0);
    let half_y = (r_sq - dx * dx).max(0.0).sqrt();
    let iy_min = (cy - half_y).max(rect.y_min);
    let iy_max = (cy + half_y).min(rect.y_max);

    // Corner intersection points
    let left_half = (r_sq - (ix_min - cx).powi(2)).max(0.0).sqrt();
    let left_y1 = (cy - left_half).max(iy_min);
    let left_y2 = (cy + left_half).min(iy_max);

    let right_half = (r_sq - (ix_max - cx).powi(2)).max(0.0).sqrt();
    let right_y1 = (cy - right_half).max(iy_min);
    let right_y2 = (cy + right_half).min(iy_max);

    let bottom_half = (r_sq - (iy_min - cy).powi(2)).max(0.0).sqrt();
    let bottom_x1 = (cx - bottom_half).max(ix_min);
    let bottom_x2 = (cx + bottom_half).min(ix_max);

    let top_half = (r_sq - (iy_max - cy).powi(2)).max(0.0).sqrt();
    let top_x1 = (cx - top_half).max(ix_min);
    let top_x2 = (cx + top_half).min(ix_max);

    let rect_area = (ix_max - ix_min) * (iy_max - iy_min);

    rect_area
        - corner_correction(bottom_x1 - ix_min, left_y1 - iy_min, r)
        - corner_correction(ix_max - bottom_x2, right_y1 - iy_min, r)
        - corner_correction(top_x1 - ix_min, iy_max - left_y2, r)
        - corner_correction(ix_max - top_x2, iy_max - right_y2, r)
}

/// Upper bound on the area between circular arc and rectangle corner.
/// The exact circular segment area is r^2 * asin(h/r) - h * sqrt(r^2 - h^2),
/// where h is half the corner diagonal. Expanding asin(x) ≈ x + x³/6 and
/// √(1-x²) ≈ 1 - x²/2 with x = h/r, the linear terms cancel and the
/// cubic terms combine to segment ≈ (2/3)*h^3/r = d^3/(12*r).
/// Dropping higher-order terms underestimates the segment
/// and overestimates this correction —
/// giving a lower bound on the circle-rect intersection area.
#[inline]
fn corner_correction(base: f32, height: f32, r: f32) -> f32 {
    if base <= 0.0 || height <= 0.0 {
        return 0.0;
    }
    let d_sq = base * base + height * height;
    let d = d_sq.sqrt();
    0.5 * base * height - d * d_sq / (12.0 * r)
}

impl QTQueryable for Circle {
    fn guarantees_collision(&self, bbox: &Rect, haz_presence_area: f32) -> bool {
        let remaining_area = bbox.area() - haz_presence_area;
        // Remaining area is too large for test to be effecient (> 0.5 * bbox.area())
        if haz_presence_area <= remaining_area {
            return false;
        }

        // Early exit: if max possible presence area can't trigger, bail
        if self.area() <= remaining_area {
            return false;
        }

        circle_rect_intersection_area(self.center.0, self.center.1, self.radius, bbox) > remaining_area
    }
}

impl QTQueryable for Rect {
    fn guarantees_collision(&self, bbox: &Rect, haz_presence_area: f32) -> bool {
        if let Some(intersection) = Rect::intersection(*self, *bbox) {
            intersection.area() > bbox.area() - haz_presence_area
        } else {
            false
        }
    }
}

impl QTQueryable for Edge {
    fn collides_with_quadrants(&self, r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        debug_assert!(r.quadrants().iter().zip(qs.iter()).all(|(q, r)| *q == **r));
        let e_x_min = self.x_min();
        let e_x_max = self.x_max();
        let e_y_min = self.y_min();
        let e_y_max = self.y_max();

        let [mut c_q0, mut c_q1, mut c_q2, mut c_q3] = [0, 1, 2, 3].map(|idx| {
            let q = qs[idx];

            let x_no_overlap = e_x_min.max(q.x_min) > e_x_max.min(q.x_max);
            let y_no_overlap = e_y_min.max(q.y_min) > e_y_max.min(q.y_max);

            if x_no_overlap || y_no_overlap {
                // Edge is completely outside the x- or y-range of the quadrant
                Some(false)
            } else if q.collides_with(&self.start) || q.collides_with(&self.end) {
                // Edge has at least one end point in the quadrant
                Some(true)
            } else {
                // Undetermined, we need to check for intersections with the sides of the quadrants
                None
            }
        });

        // If all quadrants are already determined, we can return early
        if let (Some(c_q0), Some(c_q1), Some(c_q2), Some(c_q3)) = (c_q0, c_q1, c_q2, c_q3) {
            return [c_q0, c_q1, c_q2, c_q3];
        }

        // Otherwise, we need to check for intersections with the sides of the quadrants
        // We can exploit the fact that the quadrants have a fixed layout, and share edges.

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

        //  1    0
        //  2    3

        half_intersect(self, &left, [&mut c_q1], [&mut c_q2]);
        half_intersect(self, &right, [&mut c_q3], [&mut c_q0]);
        half_intersect(self, &top, [&mut c_q0], [&mut c_q1]);
        half_intersect(self, &bottom, [&mut c_q2], [&mut c_q3]);
        half_intersect(
            self,
            &h_bisect,
            [&mut c_q1, &mut c_q2],
            [&mut c_q0, &mut c_q3],
        );
        half_intersect(
            self,
            &v_bisect,
            [&mut c_q2, &mut c_q3],
            [&mut c_q0, &mut c_q1],
        );

        let [c_q0, c_q1, c_q2, c_q3] = [c_q0, c_q1, c_q2, c_q3].map(|c| c.unwrap_or(false));
        debug_assert!(
            {
                // make sure all quadrants which are colliding according to the individual collision check are at least
                // also caught by the quadrant collision check
                qs.map(|q| self.collides_with(q))
                    .iter()
                    .zip([c_q0, c_q1, c_q2, c_q3].iter())
                    .all(|(&i_c, &q_c)| !i_c || q_c)
            },
            "{:?}, {:?}, {:?}, {:?}, {:?}",
            self,
            r,
            qs,
            [c_q0, c_q1, c_q2, c_q3],
            qs.map(|q| self.collides_with(q))
        );

        [c_q0, c_q1, c_q2, c_q3]
    }
}

impl QTQueryable for SPolygon {
    fn guarantees_collision(&self, bbox: &Rect, haz_presence_area: f32) -> bool {
        if let Some(surrogate) = &self.surrogate {
            surrogate
                .ff_poles()
                .iter()
                .any(|pole| pole.guarantees_collision(bbox, haz_presence_area))
        } else {
            false
        }
    }
}

/// If e1 intersects with e2 in the first half of e2, it sets all bools in `fst_qs` to true,
/// if e1 intersects with e2 the second half of e2, it sets all bools in `sec_qs` to true.
fn half_intersect<const N: usize>(
    e1: &Edge,
    e2: &Edge,
    fst_qs: [&mut Option<bool>; N],
    sec_qs: [&mut Option<bool>; N],
) {
    if fst_qs.iter().chain(sec_qs.iter()).any(|t| t.is_none())
        && let Some((_, e2_col_loc)) = edge_intersection_half(e1, e2)
    {
        match e2_col_loc {
            CollisionHalf::FirstHalf => {
                for c in fst_qs {
                    *c = Some(true);
                }
            }
            CollisionHalf::Halfway => {
                for c in fst_qs {
                    *c = Some(true);
                }
                for c in sec_qs {
                    *c = Some(true);
                }
            }
            CollisionHalf::SecondHalf => {
                for c in sec_qs {
                    *c = Some(true);
                }
            }
        }
    }
}

#[inline(always)]
// Similar to `edge_intersection`, but in case of an intersection, it returns in which half for both edge the intersection occurs.
fn edge_intersection_half(e1: &Edge, e2: &Edge) -> Option<(CollisionHalf, CollisionHalf)> {
    let Point(x1, y1) = e1.start;
    let Point(x2, y2) = e1.end;
    let Point(x3, y3) = e2.start;
    let Point(x4, y4) = e2.end;

    //based on: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment
    let t_nom = (x2 - x4) * (y4 - y3) - (y2 - y4) * (x4 - x3);
    let t_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
    let u_nom = (x2 - x4) * (y2 - y1) - (y2 - y4) * (x2 - x1);
    let u_denom = (x2 - x1) * (y4 - y3) - (y2 - y1) * (x4 - x3);
    if t_denom == 0.0 || u_denom == 0.0 {
        //parallel edges
        return None;
    }

    let t = t_nom / t_denom; //refers to the position along e1
    let u = u_nom / u_denom; //refers to the position along e2
    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
        let e1_loc = match t.partial_cmp(&0.5).unwrap() {
            Ordering::Greater => CollisionHalf::FirstHalf,
            Ordering::Less => CollisionHalf::SecondHalf,
            Ordering::Equal => CollisionHalf::Halfway,
        };
        let e2_loc = match u.partial_cmp(&0.5).unwrap() {
            Ordering::Greater => CollisionHalf::FirstHalf,
            Ordering::Less => CollisionHalf::SecondHalf,
            Ordering::Equal => CollisionHalf::Halfway,
        };

        return Some((e1_loc, e2_loc));
    }
    None
}

pub enum CollisionHalf {
    FirstHalf,
    Halfway,
    SecondHalf,
}
