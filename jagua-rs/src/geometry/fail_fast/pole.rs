use std::collections::VecDeque;
use std::f32::consts::PI;

use crate::geometry::geo_traits::{CollidesWith, DistanceTo};
use crate::geometry::primitives::Circle;
use crate::geometry::primitives::Rect;
use crate::geometry::primitives::SPolygon;

use anyhow::{Result, anyhow};

/// Calculate net perimeter contribution of a new pole.
/// This is the new pole's perimeter minus arcs lost to overlap (both directions).
fn net_pole_perimeter(new_pole: &Circle, existing_poles: &[Circle]) -> f32 {
    let full_perimeter = new_pole.perimeter();
    
    let lost: f32 = existing_poles.iter().map(|other| {
        new_pole.arc_inside(other) + other.arc_inside(new_pole)
    }).sum();
    
    (full_perimeter - lost).max(0.0)
}

/// Calculate net area contribution of a new pole (for coverage stopping condition).
fn net_pole_area(new_pole: &Circle, poles: &[Circle]) -> f32 {
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

/// Computes the *pole* - the circle that maximizes net perimeter contribution while staying inside `shape`.
// Poles may overlap; the perimeter-based objective naturally balances size vs. overlap to maximize effective collision surface.
/// Closely related to [Pole of Inaccessibility (PoI)](https://en.wikipedia.org/wiki/Pole_of_inaccessibility),
/// and inspired by Mapbox's [`polylabel`](https://github.com/mapbox/polylabel) algorithm.
pub fn compute_pole(shape: &SPolygon, poles: &[Circle]) -> Result<Circle> {
    let square_bbox = shape.bbox.inflate_to_square();
    let root = POINode::new(square_bbox, MAX_POI_TREE_DEPTH, shape, poles);
    let mut queue = VecDeque::from([root]);
    let mut best: Option<(Circle, f32)> = None;
    let score = |pair: &Option<(Circle, f32)>| pair.as_ref().map_or(0.0, |(_, s)| *s);

    while let Some(node) = queue.pop_front() {
        //check if better than current best
        if node.perimeter > score(&best) && node.distance > 0.0 {
            best = Some((Circle::try_new(node.bbox.centroid(), node.distance).unwrap(), node.perimeter));
        }

        //see if worth it to split
        if node.perimeter_upperbound() > score(&best)
            && let Some(children) = node.split(shape, poles)
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
    pub perimeter: f32,
}

impl POINode {
    fn new(bbox: Rect, level: usize, poly: &SPolygon, poles: &[Circle]) -> Self {
        let radius = bbox.diameter() / 2.0;
        let centroid = bbox.centroid();
        let centroid_inside = poly.collides_with(&centroid);

        let distance = {
            let distance_to_border = poly.edge_iter()
                .map(|e| e.distance_to(&centroid))
                .fold(f32::MAX, |acc, d| acc.min(d));

            //if the centroid is outside, distance is counted negative
            match centroid_inside {
                true => distance_to_border,
                false => -distance_to_border,
            }
        };

        let perimeter = if distance > 0.0 {
            let candidate = Circle { center: centroid, radius: distance };
            net_pole_perimeter(&candidate, poles)
        } else {
            0.0  // Outside polygon
        };

        Self {
            bbox,
            level,
            radius,
            distance,
            perimeter,
        }
    }

    fn split(&self, poly: &SPolygon, poles: &[Circle]) -> Option<[POINode; 4]> {
        match self.level {
            0 => None,
            _ => Some(
                self.bbox
                    .quadrants()
                    .map(|qd| POINode::new(qd, self.level - 1, poly, poles)),
            ),
        }
    }

    fn perimeter_upperbound(&self) -> f32 {
        let max_radius = (self.radius + self.distance).max(0.0);
        2.0 * PI * max_radius
    }
}
