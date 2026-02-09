use crate::collision_detection::quadtree::qt_traits::QTQueryable;
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
    pub fn from_entire_shape(shape: &SPolygon, presence_area: f32) -> Self {
        let edges = shape.edge_iter().collect();
        let ff_bbox = shape.bbox;
        let points = shape.vertices.clone();
        Self { edges, ff_bbox, points, presence_area }
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
