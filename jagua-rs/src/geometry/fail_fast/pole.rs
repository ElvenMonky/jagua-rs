use std::collections::VecDeque;
use std::f32::consts::PI;

use crate::geometry::geo_traits::{CollidesWith, DistanceTo};
use crate::geometry::primitives::Circle;
use crate::geometry::primitives::Edge;
use crate::geometry::primitives::Point;
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::SPolygon;

use anyhow::{Result, anyhow};

/// A point on a circle's boundary, stored with its angle from the circle's center
#[derive(Clone, Copy, Debug)]
struct CircleBoundaryPoint {
    point: Point,
    angle: f32,  // angle from circle center, in [-PI, PI]
}

/// Precomputed intersection data for a pole (circle)
struct PoleIntersections {
    center: Point,
    radius: f32,
    /// Boundary points sorted by angle
    boundary_points: Vec<CircleBoundaryPoint>,
}

impl PoleIntersections {
    /// Compute all intersection points of this pole with polygon edges and other poles
    fn new(pole: &Circle, poly: &SPolygon, other_poles: &[Circle]) -> Self {
        let center = pole.center;
        let radius = pole.radius;
        let mut boundary_points = Vec::new();

        // Find intersections with polygon edges
        for edge in poly.edge_iter() {
            if let Some(pt) = circle_edge_intersection(center, radius, &edge) {
                let angle = (pt.1 - center.1).atan2(pt.0 - center.0);
                boundary_points.push(CircleBoundaryPoint { point: pt, angle });
            }
        }

        // Find intersections with other poles
        for other in other_poles {
            if other.center == center && other.radius == radius {
                continue; // skip self
            }
            if let Some(pts) = circle_circle_intersections(center, radius, other.center, other.radius) {
                for pt in pts {
                    let angle = (pt.1 - center.1).atan2(pt.0 - center.0);
                    boundary_points.push(CircleBoundaryPoint { point: pt, angle });
                }
            }
        }

        // Sort by angle
        boundary_points.sort_by(|a, b| a.angle.partial_cmp(&b.angle).unwrap());

        Self {
            center,
            radius,
            boundary_points,
        }
    }

    /// For a candidate point, compute the constraint from this pole.
    /// Returns 0.0 if point is on the wrong side of the chord (blocked).
    /// Returns distance to chord if point is on the valid side.
    fn compute_constraint(&self, candidate: Point) -> f32 {
        if self.boundary_points.len() < 2 {
            // Not enough intersections to form a chord - fall back to separation distance
            let dist = self.center.distance_to(&candidate);
            return (dist - self.radius).max(0.0);
        }

        // Compute angle from pole center to candidate
        let angle = (candidate.1 - self.center.1).atan2(candidate.0 - self.center.0);

        // Find the two boundary points whose sector contains this angle
        let (p1, p2) = self.find_sector_endpoints(angle);

        // The chord is the line segment from p1 to p2
        let chord = match Edge::try_new(p1, p2) {
            Ok(e) => e,
            Err(_) => {
                // Degenerate chord - fall back to separation distance
                let dist = self.center.distance_to(&candidate);
                return (dist - self.radius).max(0.0);
            }
        };

        // Check which side of chord the candidate is on
        let dist_center_to_chord = chord.distance_to(&self.center);
        let dist_candidate_to_center = self.center.distance_to(&candidate);

        if dist_center_to_chord > dist_candidate_to_center {
            // Candidate is on the same side as center (inside the arc) - blocked
            0.0
        } else {
            // Candidate is on the far side of chord - use chord distance as constraint
            chord.distance_to(&candidate)
        }
    }

    /// Find the two boundary points whose angular sector contains the given angle
    fn find_sector_endpoints(&self, angle: f32) -> (Point, Point) {
        let n = self.boundary_points.len();
        
        // Binary search for the first point with angle > given angle
        let idx = self.boundary_points
            .binary_search_by(|bp| bp.angle.partial_cmp(&angle).unwrap())
            .unwrap_or_else(|i| i);

        // The sector is between boundary_points[idx-1] and boundary_points[idx]
        // Handle wraparound
        let idx_before = if idx == 0 { n - 1 } else { idx - 1 };
        let idx_after = idx % n;

        (self.boundary_points[idx_before].point, self.boundary_points[idx_after].point)
    }
}

/// Compute tangent point of a circle with an edge (line segment), if it exists.
/// The circle is inscribed, so it touches the edge at most at one point (tangent).
fn circle_edge_intersection(center: Point, radius: f32, edge: &Edge) -> Option<Point> {
    let closest = edge.closest_point_on_edge(&center);
    let dist_sq = (closest.0 - center.0).powi(2) + (closest.1 - center.1).powi(2);
    let radius_sq = radius * radius;

    // Check if distance from center to closest point equals radius (tangent)
    // Allow small tolerance for floating point
    if (dist_sq - radius_sq).abs() < radius_sq * 0.0001 {
        Some(closest)
    } else {
        None
    }
}

/// Compute intersection points of two circles.
/// Returns None if circles don't intersect, or Some with one or two points.
fn circle_circle_intersections(c1: Point, r1: f32, c2: Point, r2: f32) -> Option<[Point; 2]> {
    let dx = c2.0 - c1.0;
    let dy = c2.1 - c1.1;
    let d = (dx * dx + dy * dy).sqrt();

    // Check if circles intersect
    if d > r1 + r2 || d < (r1 - r2).abs() || d < 1e-10 {
        return None;
    }

    // Distance from c1 to the line through intersection points
    let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let h_sq = r1 * r1 - a * a;
    
    if h_sq < 0.0 {
        return None;
    }
    
    let h = h_sq.sqrt();

    // Point on line between centers at distance a from c1
    let px = c1.0 + a * dx / d;
    let py = c1.1 + a * dy / d;

    // Perpendicular direction
    let perp_x = -dy / d;
    let perp_y = dx / d;

    // Two intersection points (may be the same if h ≈ 0)
    Some([
        Point(px + h * perp_x, py + h * perp_y),
        Point(px - h * perp_x, py - h * perp_y),
    ])
}

/// Calculate net area contribution of a new pole (its area minus overlap with existing poles).
/// Note: This slightly underestimates net area when the new pole covers intersections
/// of existing poles, but is sufficient for the coverage stopping condition.
pub fn net_pole_area(new_pole: &Circle, poles: &[Circle]) -> f32 {
    let total_area = new_pole.area();
    let overlap: f32 = poles
        .iter()
        .map(|p| new_pole.intersection_area(p))
        .sum();
    (total_area - overlap).max(0.0)
}

///Generates a set of 'poles' for a shape according to specified coverage limits.
///See [`compute_pole`] for details on what a 'pole' is.
pub fn generate_surrogate_poles(
    shape: &SPolygon,
    n_pole_limits: &[(usize, f32)],
) -> Result<Vec<Circle>> {
    let mut all_poles = vec![shape.poi];
    let mut total_net_area = shape.poi.area();

    //Generate the poles until one of the pole number / coverage limits is reached
    loop {
        let next = compute_pole(shape, &all_poles)?;

        total_net_area += net_pole_area(&next, &all_poles);
        all_poles.push(next);

        let current_coverage = total_net_area / shape.area;

        //check if any limit in the number of poles is reached at this coverage
        let active_pole_limit = n_pole_limits
            .iter()
            .filter(|(_, coverage_threshold)| current_coverage > *coverage_threshold)
            .min_by_key(|(n_poles, _)| *n_poles)
            .map(|(n_poles, _)| n_poles);

        if let Some(active_pole_limit) = active_pole_limit
            && all_poles.len() >= *active_pole_limit
        {
            //stop generating if we are above the limit
            break;
        }
        assert!(
            all_poles.len() < 1000,
            "More than 1000 poles were generated, please check the SPSurrogateConfig"
        )
    }
    Ok(all_poles)
}

/// Computes the *pole* - the largest circle which is both inside of `shape` 
/// while respecting chord constraints from existing poles (allowing controlled overlap).
/// Closely related to [Pole of Inaccessibility (PoI)](https://en.wikipedia.org/wiki/Pole_of_inaccessibility),
/// and inspired by Mapbox's [`polylabel`](https://github.com/mapbox/polylabel) algorithm.
pub fn compute_pole(shape: &SPolygon, poles: &[Circle]) -> Result<Circle> {
    // Precompute intersection data for all existing poles
    let pole_intersections: Vec<PoleIntersections> = poles
        .iter()
        .map(|p| PoleIntersections::new(p, shape, poles))
        .collect();

    let square_bbox = shape.bbox.inflate_to_square();
    let root = POINode::new(square_bbox, MAX_POI_TREE_DEPTH, shape, poles, &pole_intersections);
    let mut queue = VecDeque::from([root]);
    let mut best: Option<(Circle, f32)> = None;
    let area = |pair: &Option<(Circle, f32)>| pair.as_ref().map_or(0.0, |(_, a)| *a);

    while let Some(node) = queue.pop_front() {
        //check if better than current best
        if node.area > area(&best) {
            best = Some((Circle::try_new(node.bbox.centroid(), node.distance).unwrap(), node.area));
        }

        //see if worth it to split
        if node.area_upperbound() > area(&best)
            && let Some(children) = node.split(shape, poles, &pole_intersections)
        {
            queue.extend(children);
        }
    }

    best.map(|(c, _)| c).ok_or(anyhow!(
        "no pole found with {} levels of recursion. Please check the input shape: {:?}",
        MAX_POI_TREE_DEPTH,
        &shape.vertices
    ))
}

const MAX_POI_TREE_DEPTH: usize = 10;

struct POINode {
    pub level: usize,
    pub bbox: Rect,
    pub radius: f32,
    pub distance: f32,
    pub area: f32,
}

impl POINode {
    fn new(bbox: Rect, level: usize, poly: &SPolygon, poles: &[Circle], pole_intersections: &[PoleIntersections]) -> Self {
        let radius = bbox.diameter() / 2.0;

        let centroid_inside = poly.collides_with(&bbox.centroid());

        let distance = {
            let distance_to_edges = poly.edge_iter().map(|e| e.distance_to(&bbox.centroid()));

            let distance_to_poles = pole_intersections
                .iter()
                .map(|pi| pi.compute_constraint(bbox.centroid()));

            let distance_to_border = distance_to_edges
                .chain(distance_to_poles)
                .fold(f32::MAX, |acc, d| acc.min(d));

            //if the centroid is outside, distance is counted negative
            match centroid_inside {
                true => distance_to_border,
                false => -distance_to_border,
            }
        };

        let area = if distance > 0.0 {
            let candidate = Circle { center: bbox.centroid(), radius: distance };
            net_pole_area(&candidate, poles)
        } else {
            0.0  // Outside polygon or blocked
        };

        Self {
            bbox,
            level,
            radius,
            distance,
            area,
        }
    }

    fn split(&self, poly: &SPolygon, poles: &[Circle], pole_intersections: &[PoleIntersections]) -> Option<[POINode; 4]> {
        match self.level {
            0 => None,
            _ => Some(
                self.bbox
                    .quadrants()
                    .map(|qd| POINode::new(qd, self.level - 1, poly, poles, pole_intersections)),
            ),
        }
    }

    fn area_upperbound(&self) -> f32 {
        let max_radius = (self.radius + self.distance).max(0.0);
        PI * max_radius * max_radius
    }
}
