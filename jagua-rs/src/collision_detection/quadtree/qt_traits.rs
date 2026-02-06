use crate::geometry::geo_traits::{CollidesWith, DistanceTo};
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::{Circle, Edge, Point, SPolygon};
use std::cmp::Ordering;
use std::f32::consts::PI;

// Quadrant bit flags
const Q0: u8 = 0b0001;
const Q1: u8 = 0b0010;
const Q2: u8 = 0b0100;
const Q3: u8 = 0b1000;
const ALL: u8 = 0b1111;

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

    /// Fast positive precheck against hazard surrogate poles.
    /// Returns true if the entity definitely collides with at least one pole.
    /// May return false negatives, never false positives.
    /// Default: no precheck (returns false).
    fn collides_with_hazard_poles(&self, _poles: &[Circle]) -> bool {
        false
    }
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

        // Count corners inside circle
        let r_sq = self.radius * self.radius;
        let corners = bbox.corners();
        let n_inside = corners.iter()
            .filter(|c| self.center.sq_distance_to(c) <= r_sq)
            .count();

        // Circle contains triangle built on 3 bbox corners, which is at least 0.5 area of the bbox
        if n_inside >= 3 {
            return true;
        }

        // radius of inscribed part of circle
        let (r, x, y) = (self.radius, self.center.0, self.center.1);
        let x0 = 0.5 * (bbox.x_max.min(x + r) + bbox.x_min.max(x - r));
        let y0 = 0.5 * (bbox.y_max.min(y + r) + bbox.y_min.max(y - r));
        let d = r - ((x - x0).powi(2) + (y - y0).powi(2)).sqrt();
        let r0 = d.min(x0 - bbox.x_min).min(y0 - bbox.y_min).min(bbox.x_max - x0).min(bbox.y_max - y0);

        PI * r0 * r0 > remaining_area
    }

    fn collides_with_hazard_poles(&self, poles: &[Circle]) -> bool {
        poles.iter().any(|pole| pole.collides_with(self))
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

    fn collides_with_hazard_poles(&self, poles: &[Circle]) -> bool {
        poles.iter().any(|pole| pole.collides_with(self))
    }
}

impl QTQueryable for Edge {
    fn collides_with_quadrants(&self, r: &Rect, qs: [&Rect; 4]) -> [bool; 4] {
        debug_assert!(r.quadrants().iter().zip(qs.iter()).all(|(q, r)| *q == **r));
        
        let mut determined: u8 = 0;
        let mut collides: u8 = 0;

        let e_x_min = self.x_min();
        let e_x_max = self.x_max();
        let e_y_min = self.y_min();
        let e_y_max = self.y_max();

        // Check bbox overlap and endpoint containment
        for (i, q) in qs.iter().enumerate() {
            let bit = 1 << i;

            let x_no_overlap = e_x_min.max(q.x_min) > e_x_max.min(q.x_max);
            let y_no_overlap = e_y_min.max(q.y_min) > e_y_max.min(q.y_max);

            if x_no_overlap || y_no_overlap {
                // Edge is completely outside the x- or y-range of the quadrant
                determined |= bit;
            } else if q.collides_with(&self.start) || q.collides_with(&self.end) {
                // Edge has at least one end point in the quadrant
                determined |= bit;
                collides |= bit;
            }
        }

        // If all quadrants are already determined, we can return early
        if determined == ALL {
            return bits_to_array(collides);
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

        half_intersect_bits(self, &left, Q1, Q2, &mut determined, &mut collides);
        half_intersect_bits(self, &right, Q3, Q0, &mut determined, &mut collides);
        half_intersect_bits(self, &top, Q0, Q1, &mut determined, &mut collides);
        half_intersect_bits(self, &bottom, Q2, Q3, &mut determined, &mut collides);
        half_intersect_bits(self, &h_bisect, Q1 | Q2, Q0 | Q3, &mut determined, &mut collides);
        half_intersect_bits(self, &v_bisect, Q2 | Q3, Q0 | Q1, &mut determined, &mut collides);

        let result = bits_to_array(collides);
        debug_assert!(
            {
                // make sure all quadrants which are colliding according to the individual collision check are at least
                // also caught by the quadrant collision check
                qs.map(|q| self.collides_with(q))
                    .iter()
                    .zip(result.iter())
                    .all(|(&i_c, &q_c)| !i_c || q_c)
            },
            "{:?}, {:?}, {:?}, {:?}, {:?}",
            self,
            r,
            qs,
            result,
            qs.map(|q| self.collides_with(q))
        );

        result
    }

    fn collides_with_hazard_poles(&self, poles: &[Circle]) -> bool {
        poles.iter().any(|pole| pole.collides_with(self))
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

#[inline(always)]
fn bits_to_array(bits: u8) -> [bool; 4] {
    [
        bits & Q0 != 0,
        bits & Q1 != 0,
        bits & Q2 != 0,
        bits & Q3 != 0,
    ]
}

/// If e1 intersects with e2 in the first half of e2, sets fst_mask bits in collides.
/// If e1 intersects with e2 in the second half of e2, sets sec_mask bits in collides.
/// Updates determined accordingly.
#[inline(always)]
fn half_intersect_bits(
    e1: &Edge,
    e2: &Edge,
    fst_mask: u8,
    sec_mask: u8,
    determined: &mut u8,
    collides: &mut u8,
) {
    let relevant = (fst_mask | sec_mask) & !*determined;
    if relevant != 0 {
        if let Some((_, e2_col_loc)) = edge_intersection_half(e1, e2) {
            let hit_mask = match e2_col_loc {
                CollisionHalf::FirstHalf => fst_mask,
                CollisionHalf::SecondHalf => sec_mask,
                CollisionHalf::Halfway => fst_mask | sec_mask,
            };
            *determined |= hit_mask;
            *collides |= hit_mask;
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
