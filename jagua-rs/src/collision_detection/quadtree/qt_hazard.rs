use crate::collision_detection::hazards::HazardEntity;
use crate::collision_detection::hazards::{HazKey, Hazard};
use crate::collision_detection::quadtree::qt_partial_hazard::QTHazPartial;
use crate::geometry::geo_enums::{GeoRelation};
use crate::geometry::primitives::{Point, Rect};
use crate::util::assertions;
use slotmap::SlotMap;
use std::array;

/// Representation of a [`Hazard`] in a [`QTNode`](crate::collision_detection::quadtree::QTNode)
#[derive(Clone, Debug)]
pub struct QTHazard {
    /// The bounding box of the quadtree node
    pub qt_bbox: Rect,
    /// The key of the hazard in the hazard map in [`CDEngine`](crate::collision_detection::cd_engine::CDEngine)
    pub hkey: HazKey,
    /// Entity inducing the hazard
    pub entity: HazardEntity,
    /// How the hazard is present in the node
    pub presence: QTHazPresence,
}

/// Presence of a [`Hazard`] in a [`QTNode`](crate::collision_detection::quadtree::QTNode)
#[derive(Clone, Debug)]
pub enum QTHazPresence {
    /// The hazard is entirely absent from the node
    None,
    /// The hazard is only partially present in the node
    Partial(QTHazPartial),
    /// The hazard is present in the entire node.
    Entire,
}
impl QTHazard {
    /// Converts a [`Hazard`] into a [`QTHazard`], assuming it is for the root of the quadtree.
    pub fn from_root(qt_root_bbox: Rect, haz: &Hazard, hkey: HazKey) -> Self {
        Self {
            qt_bbox: qt_root_bbox,
            hkey,
            entity: haz.entity,
            presence: QTHazPresence::Partial(QTHazPartial::from_entire_shape(&haz.shape, haz.shape.area)),
        }
    }

    /// Returns the resulting QTHazards after constricting to the provided quadrants.
    /// The quadrants should be ordered according to the [Cartesian system](https://en.wikipedia.org/wiki/Quadrant_(plane_geometry))
    /// and should all be inside the bounds from which `self` was created.
    pub fn constrict(&self, quadrants: [Rect; 4], haz_map: &SlotMap<HazKey, Hazard>) -> [Self; 4] {
        debug_assert!(
            quadrants
                .iter()
                .all(|q| self.qt_bbox.relation_to(*q) == GeoRelation::Surrounding)
        );
        debug_assert!(assertions::quadrants_have_valid_layout(&quadrants));

        match &self.presence {
            QTHazPresence::None => unreachable!("Hazard presence cannot be None in a QTHazard"),
            QTHazPresence::Entire => array::from_fn(|_| self.clone()), // The hazard is entirely present in all quadrants
            QTHazPresence::Partial(partial_haz) => {
                //If the hazard is partially present, we need to check which type of presence each quadrant has

                let haz_shape = &haz_map[self.hkey].shape;

                //Check if one of the quadrants entirely contains the hazard
                let enclosed_hazard_quadrant = quadrants
                    .iter()
                    .map(|q| haz_shape.bbox.relation_to(*q))
                    .position(|r| r == GeoRelation::Enclosed);

                if let Some(quad_index) = enclosed_hazard_quadrant {
                    //The hazard is entirely enclosed within one quadrant,
                    //For this quadrant the QTHazard is equivalent to the original hazard, the rest are None
                    return array::from_fn(|i| {
                        let presence = if i == quad_index {
                            self.presence.clone()
                        } else {
                            QTHazPresence::None
                        };
                        Self {
                            qt_bbox: quadrants[i],
                            presence,
                            hkey: self.hkey,
                            entity: self.entity,
                        }
                    });
                }
                //The hazard is active in multiple quadrants

                let center = self.qt_bbox.centroid();
                let cx = center.0;
                let cy = center.1;

                //  Q1 | Q0
                //  ---+---
                //  Q2 | Q3
                //
                // Quadrant index from point position.
                // Points exactly on a midline stay in `prev` quadrant,
                // preferring adjacent over opposite.
                let quadrant_of = |p: Point, prev: usize| -> usize {
                    let on_v = p.0 == cx;
                    let on_h = p.1 == cy;
                    if on_v && on_h {
                        return prev;
                    }
                    let right = if on_v { prev == 0 || prev == 3 } else { p.0 > cx };
                    let upper = if on_h { prev == 0 || prev == 1 } else { p.1 > cy };
                    match (right, upper) {
                        (true, true) => 0,
                        (false, true) => 1,
                        (false, false) => 2,
                        (true, false) => 3,
                    }
                };

                let cross_v = |p1: Point, p2: Point| -> Point {
                    let t = (cx - p1.0) / (p2.0 - p1.0);
                    Point(cx, p1.1 + t * (p2.1 - p1.1))
                };
                let cross_h = |p1: Point, p2: Point| -> Point {
                    let t = (cy - p1.1) / (p2.1 - p1.1);
                    Point(p1.0 + t * (p2.0 - p1.0), cy)
                };

                // Working storage for the 4 quadrant partial hazards
                let mut q_haz: [QTHazPartial; 4] = array::from_fn(|_| QTHazPartial {
                    edges: Vec::new(),
                    ff_bbox: partial_haz.ff_bbox,
                    points: Vec::new(),
                    presence_area: 0.0,
                });

                // === Single pass over edges ===
                for &edge in &partial_haz.edges {
                    let s_on_v = edge.start.0 == cx;
                    let s_on_h = edge.start.1 == cy;
                    let e_on_v = edge.end.0 == cx;
                    let e_on_h = edge.end.1 == cy;
                    let s_on_midline = s_on_v || s_on_h;
                    let e_on_midline = e_on_v || e_on_h;

                    let (sq, eq) = if !s_on_midline && !e_on_midline {
                        // Neither on midline — standard case
                        let sq = quadrant_of(edge.start, 0);
                        let eq = quadrant_of(edge.end, sq);
                        (sq, eq)
                    } else if !e_on_midline {
                        // Only start on midline — anchor from end
                        let eq = quadrant_of(edge.end, 0);
                        let sq = quadrant_of(edge.start, eq);
                        (sq, eq)
                    } else if !s_on_midline {
                        // Only end on midline — anchor from start
                        let sq = quadrant_of(edge.start, 0);
                        let eq = quadrant_of(edge.end, sq);
                        (sq, eq)
                    } else if s_on_v && e_on_v {
                        // Both on vertical midline
                        if edge.start.1 > edge.end.1 {
                            // Going down → interior right (Q0/Q3)
                            let sq = if edge.start.1 > cy { 0 } else { 3 };
                            let eq = if edge.end.1 < cy { 3 } else { 0 };
                            (sq, eq)
                        } else {
                            // Going up → interior left (Q1/Q2)
                            let sq = if edge.start.1 < cy { 2 } else { 1 };
                            let eq = if edge.end.1 > cy { 1 } else { 2 };
                            (sq, eq)
                        }
                    } else if s_on_h && e_on_h {
                        // Both on horizontal midline
                        if edge.start.0 > edge.end.0 {
                            // Going left → interior below (Q2/Q3)
                            let sq = if edge.start.0 > cx { 3 } else { 2 };
                            let eq = if edge.end.0 < cx { 2 } else { 3 };
                            (sq, eq)
                        } else {
                            // Going right → interior above (Q0/Q1)
                            let sq = if edge.start.0 < cx { 1 } else { 0 };
                            let eq = if edge.end.0 > cx { 0 } else { 1 };
                            (sq, eq)
                        }
                    } else {
                        // Ends on different midlines — shared corner quadrant
                        let mid_x = if s_on_v { edge.end.0 } else { edge.start.0 };
                        let mid_y = if s_on_h { edge.end.1 } else { edge.start.1 };
                        let q = quadrant_of(Point(mid_x, mid_y), 0);
                        (q, q)
                    };

                    q_haz[sq].edges.push(edge);

                    let rel = sq ^ eq;
                    if rel == 2 {
                        let dx = edge.end.0 - edge.start.0;
                        let dy = edge.end.1 - edge.start.1;
                        let cross = dx * (cy - edge.start.1) - dy * (cx - edge.start.0);
                        if cross > 0.0 {
                            q_haz[(sq + 1) & 3].edges.push(edge);
                        } else if cross < 0.0 {
                            q_haz[(sq + 3) & 3].edges.push(edge);
                        }
                    }
                    if rel != 0 {
                        q_haz[eq].edges.push(edge);
                    }
                }

                // === Single pass over vertices for points + area ===
                let n = partial_haz.points.len();
                debug_assert!(n > 0, "Partial hazard with no points");
                let mut prev_p = partial_haz.points[n - 1];
                let mut cur_q = quadrant_of(prev_p, 0);

                macro_rules! push {
                    ($qi:expr, $p:expr) => {{
                        let qi = $qi;
                        let p = $p;
                        let h = &mut q_haz[qi];
                        if h.points.last() != Some(&p) {
                            if let Some(&last) = h.points.last() {
                                if qi != cur_q && last.0 != p.0 && last.1 != p.1 {
                                    h.presence_area += (last.1 + cy) * (last.0 - cx);
                                    h.points.push(center);
                                    h.presence_area += (cy + p.1) * (cx - p.0);
                                } else {
                                    h.presence_area += (last.1 + p.1) * (last.0 - p.0);
                                }
                            }
                            h.points.push(p);
                            cur_q = qi;
                        }
                    }};
                }

                for i in 0..n {
                    let p = partial_haz.points[i];
                    let new_q = quadrant_of(p, cur_q);

                    if new_q != cur_q {
                        let rel = cur_q ^ new_q;
                        if rel == 1 {
                            let cp = cross_v(prev_p, p);
                            push!(cur_q, cp);
                            push!(new_q, cp);
                        } else if rel == 3 {
                            let cp = cross_h(prev_p, p);
                            push!(cur_q, cp);
                            push!(new_q, cp);
                        } else {
                            let dx = p.0 - prev_p.0;
                            let dy = p.1 - prev_p.1;
                            let cross = dx * (cy - prev_p.1) - dy * (cx - prev_p.0);

                            if cross > 0.0 {
                                let mid_q = (cur_q + 1) & 3;
                                let cp1 = cross_v(prev_p, p);
                                let cp2 = cross_h(prev_p, p);
                                push!(cur_q, cp1);
                                push!(mid_q, cp1);
                                push!(mid_q, cp2);
                                push!(new_q, cp2);
                            } else if cross < 0.0 {
                                let cp1 = cross_h(prev_p, p);
                                let cp2 = cross_v(prev_p, p);
                                let mid_q = (cur_q + 3) & 3;
                                push!(cur_q, cp1);
                                push!(mid_q, cp1);
                                push!(mid_q, cp2);
                                push!(new_q, cp2);
                            } else {
                                push!(cur_q, center);
                                push!(new_q, center);
                            }
                        }
                    }

                    push!(cur_q, p);
                    prev_p = p;
                }

                // === Finalize and build results ===
                let quadrant_area = quadrants[0].area();

                array::from_fn(|i| {
                    // Close shoelace
                    let pts = &q_haz[i].points;
                    if pts.len() >= 2 {
                        let last = *pts.last().unwrap();
                        let first = pts[0];
                        q_haz[i].presence_area += (last.1 + first.1) * (last.0 - first.0);
                    }
                    q_haz[i].presence_area = (0.5 * q_haz[i].presence_area).abs();

                    let presence = if q_haz[i].presence_area > quadrant_area - 1e-9 {
                        QTHazPresence::Entire
                    } else if q_haz[i].points.is_empty() {
                        QTHazPresence::None
                    } else {
                        // Compute tight ff_bbox
                        if q_haz[i].edges.len() != partial_haz.edges.len() {
                            let (mut x_min, mut y_min, mut x_max, mut y_max) = (
                                f32::INFINITY, f32::INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY,
                            );
                            for edge in &q_haz[i].edges {
                                x_min = x_min.min(edge.start.x()).min(edge.end.x());
                                y_min = y_min.min(edge.start.y()).min(edge.end.y());
                                x_max = x_max.max(edge.start.x()).max(edge.end.x());
                                y_max = y_max.max(edge.start.y()).max(edge.end.y());
                            }
                            if x_min < x_max && y_min < y_max {
                                q_haz[i].ff_bbox = Rect { x_min, y_min, x_max, y_max };
                            }
                        }
                        QTHazPresence::Partial(std::mem::replace(
                            &mut q_haz[i],
                            QTHazPartial {
                                edges: Vec::new(),
                                ff_bbox: partial_haz.ff_bbox,
                                points: Vec::new(),
                                presence_area: 0.0,
                            },
                        ))
                    };

                    Self {
                        qt_bbox: quadrants[i],
                        presence,
                        hkey: self.hkey,
                        entity: self.entity,
                    }
                })
            }
        }
    }

    pub fn n_edges(&self) -> usize {
        match &self.presence {
            QTHazPresence::None | QTHazPresence::Entire => 0,
            QTHazPresence::Partial(partial_haz) => partial_haz.n_edges(),
        }
    }
}
