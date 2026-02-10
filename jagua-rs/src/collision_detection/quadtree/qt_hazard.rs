use crate::collision_detection::hazards::HazardEntity;
use crate::collision_detection::hazards::{HazKey, Hazard};
use crate::collision_detection::quadtree::qt_partial_hazard::QTHazPartial;
use crate::geometry::geo_enums::{GeoPosition, GeoRelation};
use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::{Point, Rect, SPolygon};
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
            presence: QTHazPresence::Partial(QTHazPartial::from_entire_shape(
                &haz.shape, haz.shape.area, haz.entity.scope(), qt_root_bbox
            )),
        }
    }

    /// Returns the resulting QTHazards after constricting to the provided quadrants.
    /// The quadrants should be ordered according to the [Cartesian system](https://en.wikipedia.org/wiki/Quadrant_(plane_geometry))
    /// and should all be inside the bounds from which `self` was created.
    pub fn constrict_old(&self, quadrants: [Rect; 4], haz_map: &SlotMap<HazKey, Hazard>) -> [Self; 4] {
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
                    array::from_fn(|i| {
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
                    })
                } else {
                    //The hazard is active in multiple quadrants

                    // First lets find the quadrants where edges of the partial hazard are colliding with the quadrants.
                    // These will also be partially present hazards.
                    let mut constricted_hazards = quadrants.map(|q| {
                        //For every quadrant, collect the edges that are colliding with it
                        let mut colliding_edges = None;
                        for edge in partial_haz.edges.iter() {
                            if q.collides_with(edge) {
                                colliding_edges.get_or_insert_with(Vec::new).push(*edge);
                            }
                        }
                        //If there are relevant edges, create a new QTHazard for this quadrant which is partially present
                        colliding_edges.map(|edges| {
                            // Compute clamped points for this quadrant
                            let mut points: Vec<Point> = Vec::new();
                            for vert in &partial_haz.points {
                                let x = vert.0.clamp(q.x_min, q.x_max);
                                let y = vert.1.clamp(q.y_min, q.y_max);
                                let p = Point(x, y);

                                let len = points.len();
                                if len >= 2 {
                                    // Collinear on vertical line - update last point's y
                                    if points[len - 1].0 == x && points[len - 2].0 == x {
                                        points[len - 1].1 = y;
                                    // Collinear on horizontal line - update last point's x
                                    } else if points[len - 1].1 == y && points[len - 2].1 == y {
                                        points[len - 1].0 = x;
                                    } else if points.last() != Some(&p) {
                                        points.push(p);
                                    }
                                } else if points.last() != Some(&p) {
                                    points.push(p);
                                }
                            }

                            // Clean up wrap-around collinearity
                            while points.len() > 1 && points.first() == points.last() {
                                points.pop();
                            }
                            if points.len() >= 3 {
                                let len = points.len();
                                if (points[0].0 == points[len - 1].0 && points[len - 1].0 == points[len - 2].0)
                                    || (points[0].1 == points[len - 1].1 && points[len - 1].1 == points[len - 2].1)
                                {
                                    points.pop();
                                }
                            }

                            let mut presence_area = SPolygon::calculate_area(&points).abs();
                            if self.entity.scope() == GeoPosition::Exterior {
                                presence_area = quadrants[0].area() - presence_area
                            };
                            QTHazard {
                                qt_bbox: q,
                                presence: QTHazPresence::Partial(QTHazPartial::from_parent(
                                    partial_haz,
                                    edges,
                                    points,
                                    presence_area,
                                )),
                                hkey: self.hkey,
                                entity: self.entity,
                            }
                        })
                    });

                    //At this point, we have resolved all quadrants that have edges colliding with them (i.e. `Partial` presence).
                    //What remain are the quadrants without any intersecting edges.
                    //These can either have the hazard `Entire` or `None` presence
                    for i in 0..4 {
                        let quadrant = quadrants[i];
                        if constricted_hazards[i].is_none() {
                            let colliding = haz_shape.collides_with(&quadrant.centroid());
                            let presence = match self.entity.scope() {
                                GeoPosition::Interior if colliding => QTHazPresence::Entire,
                                GeoPosition::Exterior if !colliding => {
                                    QTHazPresence::Entire
                                }
                                _ => QTHazPresence::None,
                            };

                            constricted_hazards[i] = Some(QTHazard {
                                qt_bbox: quadrant,
                                presence,
                                hkey: self.hkey,
                                entity: self.entity,
                            });
                        }
                    }

                    constricted_hazards
                        .map(|h| h.expect("all constricted hazards should be resolved"))
                }
            }
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
                let mut cur_q = quadrant_of(partial_haz.points[0], 0);
                cur_q = quadrant_of(prev_p, cur_q);

                macro_rules! push {
                    ($qi:expr, $p:expr, $a:expr) => {{
                        let qi = $qi;
                        let p = $p;
                        let a = $a;
                        assert!(
                            p.0.is_finite() && p.1.is_finite() &&
                            p.0 >= quadrants[qi].x_min - 1e-3 && p.0 <= quadrants[qi].x_max + 1e-3 &&
                            p.1 >= quadrants[qi].y_min - 1e-3 && p.1 <= quadrants[qi].y_max + 1e-3,
                            "{:?} OOB point {:?} pushed to quadrant {a} {cx} {cy} {qi} ({:?})\nprev_p={:?} cur_q={}\nparent points: {:?}\nparent edges: {:?}\nhaz vertices: {:?}",
                            self.qt_bbox,
                            p, quadrants[qi], prev_p, cur_q,
                            partial_haz.points, partial_haz.edges,
                            q_haz[qi].points
                        );
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
                            } else if qi != cur_q {
                                h.points.push(center);
                                h.presence_area += (cy + p.1) * (cx - p.0);
                            }
                            h.points.push(p);
                        }
                        cur_q = qi;
                    }};
                }

                for i in 0..n {
                    let p = partial_haz.points[i];
                    let new_q = quadrant_of(p, cur_q);

                    if new_q != cur_q {
                        let rel = cur_q ^ new_q;
                        if rel == 1 {
                            let t = (cx - prev_p.0) / (p.0 - prev_p.0);
                            let cp = Point(cx, prev_p.1 + t * (p.1 - prev_p.1));
                            push!(cur_q, cp, 1);
                            push!(new_q, cp, 2);
                        } else if rel == 3 {
                            let t = (cy - prev_p.1) / (p.1 - prev_p.1);
                            let cp = Point(prev_p.0 + t * (p.0 - prev_p.0), cy);
                            push!(cur_q, cp, 3);
                            push!(new_q, cp, 4);
                        } else {
                            let tv = (cx - prev_p.0) / (p.0 - prev_p.0);
                            let th = (cy - prev_p.1) / (p.1 - prev_p.1);
                            if th == tv {
                                push!(cur_q, center, 13);
                                push!(new_q, center, 14);
                            } else {
                                let cpv = Point(cx, prev_p.1 + tv * (p.1 - prev_p.1));
                                let cph = Point(prev_p.0 + th * (p.0 - prev_p.0), cy);
                                let (cp1, cp2, mid_q) = if tv < th {
                                    (cpv, cph, cur_q ^ 1)
                                } else {
                                    (cph, cpv, cur_q ^ 3)
                                };
                                push!(cur_q, cp1, 5);
                                push!(mid_q, cp1, 6);
                                push!(mid_q, cp2, 7);
                                push!(new_q, cp2, 8);
                            }
                        }
                    }

                    push!(cur_q, p, 9);
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
                    q_haz[i].presence_area = if self.entity.scope() == GeoPosition::Exterior {
                        quadrant_area - (0.5 * q_haz[i].presence_area).abs()
                    } else {
                        (0.5 * q_haz[i].presence_area).abs()
                    };

                    let presence = if q_haz[i].presence_area >= quadrant_area - 1e-4 {
                        QTHazPresence::Entire
                    } else if q_haz[i].points.is_empty() || q_haz[i].presence_area <= 1e-4 {
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

    pub fn constrict_new(&self, quadrants: [Rect; 4], haz_map: &SlotMap<HazKey, Hazard>) -> [Self; 4] {
        let result = self.constrict(quadrants, haz_map);
        
        #[cfg(debug_assertions)]
        {
            let old_result = self.constrict_old(quadrants, haz_map);
            for i in 0..4 {
                let old_p = &old_result[i].presence;
                let new_p = &result[i].presence;
                assert_eq!(old_result[i].hkey, result[i].hkey);
                assert_eq!(old_result[i].entity, result[i].entity);
                assert_eq!(old_result[i].qt_bbox, result[i].qt_bbox);
                match (old_p, new_p) {
                    (QTHazPresence::None, QTHazPresence::None) => {},
                    (QTHazPresence::Entire, QTHazPresence::Entire) => {},
                    (QTHazPresence::Partial(old), QTHazPresence::Partial(new)) => {
                        assert_eq!(old.ff_bbox, new.ff_bbox,
                            "ff_bbox mismatch in quadrant {i}");
                        assert!(old.edges.len() <= new.edges.len(),
                            "edge count mismatch in quadrant {i}: {old_p:?} vs {new_p:?} {:?} {:?}", haz_map[result[i].hkey], result[i].qt_bbox);
                        assert!(old.edges.len() + 1 >= new.edges.len(),
                            "edge count mismatch in quadrant {i}: {old_p:?} vs {new_p:?} {:?} {:?}", haz_map[result[i].hkey], result[i].qt_bbox);
                        // edges may be in different order, so compare as sets
                        // assert_eq!(old.edges, new.edges,
                        //    "edges mismatch in quadrant {i}");
                        //assert_eq!(old.points.len(), new.points.len(),
                        //    "point count mismatch in quadrant {i}: {} vs {}", old.points.len(), new.points.len());
                        // edges may be in different order, so compare as sets
                        //assert_eq!(old.points, new.points,
                        //    "points mismatch in quadrant {i}");
                        // presence_area should be close
                        assert!(old.presence_area * 0.95 < new.presence_area,
                            "area mismatch in quadrant {i}: {old_p:?} vs {new_p:?} {:?} {:?}", haz_map[result[i].hkey], result[i].qt_bbox);
                        assert!(old.presence_area * 1.05 > new.presence_area,
                            "area mismatch in quadrant {i}: {old_p:?} vs {new_p:?} {:?} {:?}", haz_map[result[i].hkey], result[i].qt_bbox);
                    },
                    _ => panic!("presence type mismatch in quadrant {i}: {old_p:?} vs {new_p:?} {:?} {:?}", haz_map[result[i].hkey], result[i].qt_bbox),
                }
            }
        }
        
        result
    }

    pub fn n_edges(&self) -> usize {
        match &self.presence {
            QTHazPresence::None | QTHazPresence::Entire => 0,
            QTHazPresence::Partial(partial_haz) => partial_haz.n_edges(),
        }
    }
}
