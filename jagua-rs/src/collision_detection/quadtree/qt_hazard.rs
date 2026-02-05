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
        let presence = haz.shape.area / qt_root_bbox.area();
        Self {
            qt_bbox: qt_root_bbox,
            hkey,
            entity: haz.entity,
            presence: QTHazPresence::Partial(QTHazPartial::from_entire_shape(&haz.shape, presence)),
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
                    let new_presence = haz_shape.area / quadrants[quad_index].area();
                    array::from_fn(|i| {
                        let presence = if i == quad_index {
                            QTHazPresence::Partial(QTHazPartial::from_parent(
                                partial_haz,
                                partial_haz.edges.clone(),
                                partial_haz.points.clone(),
                                new_presence,
                            ))
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
                    
                    let x_mid = quadrants[0].x_min;
                    let y_mid = quadrants[0].y_min;

                    // Clamping directions per quadrant: (x_high, y_high)
                    const CLAMP_DIRS: [(bool, bool); 4] = [
                        (true, true),   // Q0 (++)
                        (false, true),  // Q1 (-+)
                        (false, false), // Q2 (--)
                        (true, false),  // Q3 (+-)
                    ];

                    // First, collect edges and compute clamped points for each quadrant
                    let mut constricted_hazards = quadrants.map(|q| {
                        //For every quadrant, collect the edges that are colliding with it
                        let mut colliding_edges = None;
                        for edge in partial_haz.edges.iter() {
                            if q.collides_with(edge) {
                                colliding_edges.get_or_insert_with(Vec::new).push(*edge);
                            }
                        }
                        colliding_edges
                    });

                    // Compute clamped points for presence calculation
                    let mut q_points: [Vec<Point>; 4] = Default::default();
                    for vert in &partial_haz.points {
                        for q in 0..4 {
                            let (x_high, y_high) = CLAMP_DIRS[q];
                            let x = if x_high { vert.0.max(x_mid) } else { vert.0.min(x_mid) };
                            let y = if y_high { vert.1.max(y_mid) } else { vert.1.min(y_mid) };
                            let p = Point(x, y);

                            let len = q_points[q].len();
                            if len >= 2 {
                                let last = q_points[q][len - 1];
                                let prev = q_points[q][len - 2];
                                // Collinear on vertical line - update last point's y
                                if last.0 == x && prev.0 == x {
                                    q_points[q][len - 1].1 = y;
                                // Collinear on horizontal line - update last point's x
                                } else if last.1 == y && prev.1 == y {
                                    q_points[q][len - 1].0 = x;
                                } else if q_points[q].last() != Some(&p) {
                                    q_points[q].push(p);
                                }
                            } else if q_points[q].last() != Some(&p) {
                                q_points[q].push(p);
                            }
                        }
                    }

                    // Clean up wrap-around collinearity
                    for verts in &mut q_points {
                        while verts.len() > 1 && verts.first() == verts.last() {
                            verts.pop();
                        }
                        if verts.len() >= 3 {
                            let first = verts[0];
                            let last = verts[verts.len() - 1];
                            let prev = verts[verts.len() - 2];
                            if (first.0 == last.0 && last.0 == prev.0) || (first.1 == last.1 && last.1 == prev.1) {
                                verts.pop();
                            }
                        }
                    }

                    // Build QTHazards from edges and points
                    let mut constricted_hazards: [Option<QTHazard>; 4] = std::array::from_fn(|i| {
                        constricted_hazards[i].take().map(|edges| {
                            let presence = SPolygon::calculate_area(&q_points[i]).abs() / quadrants[i].area();
                            let points = std::mem::take(&mut q_points[i]);
                            QTHazard {
                                qt_bbox: quadrants[i],
                                presence: QTHazPresence::Partial(QTHazPartial::from_parent(
                                    partial_haz,
                                    edges,
                                    points,
                                    presence,
                                )),
                                hkey: self.hkey,
                                entity: self.entity,
                            }
                        })
                    });

                    debug_assert!(constricted_hazards.iter().filter(|h| h.is_some()).count() > 0);

                    //At this point, we have resolved all quadrants that have edges colliding with them (i.e. `Partial` presence).
                    //What remain are the quadrants without any intersecting edges.
                    //These can either have the hazard `Entire` or `None` presence
                    for i in 0..4 {
                        let quadrant = quadrants[i];
                        if constricted_hazards[i].is_none() {
                            //One important observation is that `Entire` and `None` present hazards will always be separated by a node with `Partial` presence.
                            //If a neighbor is already resolved to `Entire` or `None`, this quadrant will have the same presence.
                            //This saves quite a bit of containment checks.

                            let neighbor_presences = Rect::QUADRANT_NEIGHBOR_LAYOUT[i]
                                .map(|idx| constricted_hazards[idx].as_ref().map(|h| &h.presence));

                            let none_neighbor = neighbor_presences
                                .iter()
                                .flatten()
                                .any(|p| matches!(p, QTHazPresence::None));
                            let entire_neighbor = neighbor_presences
                                .iter()
                                .flatten()
                                .any(|p| matches!(p, QTHazPresence::Entire));

                            let presence = match (none_neighbor, entire_neighbor) {
                                (true, true) => unreachable!(
                                    "No unresolved quadrant should not have both None and Entire neighbors, this indicates a bug in the quadtree construction logic."
                                ),
                                (true, false) => QTHazPresence::None,
                                (false, true) => QTHazPresence::Entire,
                                (false, false) => {
                                    let colliding = haz_shape.collides_with(&quadrant.centroid());
                                    match self.entity.scope() {
                                        GeoPosition::Interior if colliding => QTHazPresence::Entire,
                                        GeoPosition::Exterior if !colliding => {
                                            QTHazPresence::Entire
                                        }
                                        _ => QTHazPresence::None,
                                    }
                                }
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

    pub fn n_edges(&self) -> usize {
        match &self.presence {
            QTHazPresence::None | QTHazPresence::Entire => 0,
            QTHazPresence::Partial(partial_haz) => partial_haz.n_edges(),
        }
    }
}
