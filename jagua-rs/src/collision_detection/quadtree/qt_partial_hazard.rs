use crate::collision_detection::quadtree::qt_traits::QTQueryable;
use crate::geometry::geo_enums::GeoPosition;
use crate::geometry::geo_traits::CollidesWith;
use crate::geometry::primitives::{Edge, Point, Rect, SPolygon};

/// Defines a set of edges from a hazard that is partially active in the [`QTNode`](crate::collision_detection::quadtree::QTNode).
#[derive(Clone, Debug)]
pub struct QTHazPartial {
    /// The edges that are active in the quadtree-node.
    pub edges: Vec<Edge>,
    /// A bounding box that guarantees all edges are contained within it. (used for fail fast)
    pub ff_bbox: Rect,
    /// Deduplicated outline points of clamped hazard polygon
    pub points: Vec<Point>,
    /// Lower bound on the fraction of the quadtree node's area blocked by this hazard.
    pub presence_area: f32,
}

impl QTHazPartial {
    pub fn from_entire_shape(shape: &SPolygon, presence_area: f32, scope: GeoPosition, bbox: Rect) -> Self {
        let mut edges: Vec<Edge> = shape.edge_iter().collect();
        let ff_bbox = shape.bbox;
        let mut points: Vec<Point> = Vec::new();
        for vert in shape.vertices.iter() {
            let x = vert.0.clamp(bbox.x_min, bbox.x_max);
            let y = vert.1.clamp(bbox.y_min, bbox.y_max);
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
        if scope == GeoPosition::Exterior {
            points.reverse();
            edges = edges.into_iter().map(|e| e.reverse()).collect();
        }
        Self { edges, ff_bbox, points, presence_area }
    }

    pub fn from_parent(parent: &QTHazPartial, restricted_edges: Vec<Edge>, points: Vec<Point>, presence_area: f32) -> Self {
        debug_assert!(!restricted_edges.is_empty());
        debug_assert!(restricted_edges.iter().all(|e| parent.edges.contains(e)));
        let ff_bbox = {
            //calculate a bounding box around the edges
            if parent.edges.len() == restricted_edges.len() {
                // If the edges cover the entire shape, use the shape's bounding box
                parent.ff_bbox
            } else {
                // Otherwise, calculate the bounding box from the edges
                let (mut x_min, mut y_min, mut x_max, mut y_max) = (
                    f32::INFINITY,
                    f32::INFINITY,
                    f32::NEG_INFINITY,
                    f32::NEG_INFINITY,
                );
                for edge in &restricted_edges {
                    x_min = x_min.min(edge.start.x()).min(edge.end.x());
                    y_min = y_min.min(edge.start.y()).min(edge.end.y());
                    x_max = x_max.max(edge.start.x()).max(edge.end.x());
                    y_max = y_max.max(edge.start.y()).max(edge.end.y());
                }
                if x_min < x_max && y_min < y_max {
                    Rect {
                        x_min,
                        y_min,
                        x_max,
                        y_max,
                    }
                } else {
                    // If the edges are all aligned to an axis, use the shape's bounding box
                    parent.ff_bbox
                }
            }
        };

        Self {
            edges: restricted_edges,
            ff_bbox,
            points,
            presence_area,
        }
    }

    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }
}

impl<T: QTQueryable> CollidesWith<T> for QTHazPartial {
    fn collides_with(&self, entity: &T) -> bool {
        // If the entity does not collide with the bounding box of the hazard, it cannot collide with the hazard
        entity.collides_with(&self.ff_bbox) && self.edges.iter().any(|e| entity.collides_with(e))
    }
}
