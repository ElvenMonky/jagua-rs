use crate::geometry::primitives::Circle;
use crate::geometry::primitives::Point;
use crate::geometry::geo_traits::DistanceTo;

/// Computes the smallest enclosing circle (`bounding_circle`) of a set of points.
/// Uses Welzl's algorithm, which runs in expected O(n) time.
///
/// The resulting circle fully contains all points.
pub fn smallest_enclosing_circle(points: &[Point]) -> Circle {
    assert!(!points.is_empty(), "Cannot compute smallest enclosing circle of empty point set");

    // Shuffle for expected linear time — use a deterministic permutation
    // via a simple index-reversal to avoid requiring rand dependency
    let mut indices: Vec<usize> = (0..points.len()).collect();
    // Simple deterministic shuffle: reverse pairs from outside in
    let n = indices.len();
    for i in 0..n / 2 {
        if i % 2 == 1 {
            indices.swap(i, n - 1 - i);
        }
    }

    let mut circle = circle_from_one(points[indices[0]]);

    for i in 1..indices.len() {
        let p = points[indices[i]];
        if !circle_contains(&circle, p) {
            circle = welzl_with_one(points, &indices[..i], p);
        }
    }

    circle
}

fn welzl_with_one(points: &[Point], indices: &[usize], b1: Point) -> Circle {
    let mut circle = circle_from_one(b1);

    for (j, &idx) in indices.iter().enumerate() {
        let p = points[idx];
        if !circle_contains(&circle, p) {
            circle = welzl_with_two(points, &indices[..j], b1, p);
        }
    }

    circle
}

fn welzl_with_two(points: &[Point], indices: &[usize], b1: Point, b2: Point) -> Circle {
    let mut circle = circle_from_two(b1, b2);

    for &idx in indices {
        let p = points[idx];
        if !circle_contains(&circle, p) {
            circle = circle_from_three(b1, b2, p);
        }
    }

    circle
}

fn circle_contains(c: &Circle, p: Point) -> bool {
    // Small tolerance for floating point
    c.center.sq_distance_to(&p) <= c.radius * c.radius * (1.0 + 1e-6)
}

fn circle_from_one(p: Point) -> Circle {
    Circle { center: p, radius: 0.0 }
}

fn circle_from_two(p1: Point, p2: Point) -> Circle {
    let center = Point(
        (p1.0 + p2.0) / 2.0,
        (p1.1 + p2.1) / 2.0,
    );
    let radius = center.distance_to(&p1);
    Circle { center, radius }
}

fn circle_from_three(p1: Point, p2: Point, p3: Point) -> Circle {
    // Circumscribed circle of triangle p1-p2-p3
    let ax = p1.0;
    let ay = p1.1;
    let bx = p2.0;
    let by = p2.1;
    let cx = p3.0;
    let cy = p3.1;

    let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));

    if d.abs() < 1e-10 {
        // Collinear points — fall back to the two farthest apart
        let d12 = p1.sq_distance_to(&p2);
        let d23 = p2.sq_distance_to(&p3);
        let d13 = p1.sq_distance_to(&p3);
        if d12 >= d23 && d12 >= d13 {
            return circle_from_two(p1, p2);
        } else if d23 >= d13 {
            return circle_from_two(p2, p3);
        } else {
            return circle_from_two(p1, p3);
        }
    }

    let ux = ((ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by))
        / d;
    let uy = ((ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax))
        / d;

    let center = Point(ux, uy);
    let radius = center.distance_to(&p1);
    Circle { center, radius }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_point() {
        let c = smallest_enclosing_circle(&[Point(1.0, 2.0)]);
        assert_eq!(c.center, Point(1.0, 2.0));
        assert_eq!(c.radius, 0.0);
    }

    #[test]
    fn test_two_points() {
        let c = smallest_enclosing_circle(&[Point(0.0, 0.0), Point(2.0, 0.0)]);
        assert!((c.center.0 - 1.0).abs() < 1e-5);
        assert!((c.center.1 - 0.0).abs() < 1e-5);
        assert!((c.radius - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_square() {
        let c = smallest_enclosing_circle(&[
            Point(0.0, 0.0),
            Point(1.0, 0.0),
            Point(1.0, 1.0),
            Point(0.0, 1.0),
        ]);
        assert!((c.center.0 - 0.5).abs() < 1e-4);
        assert!((c.center.1 - 0.5).abs() < 1e-4);
        let expected_r = (0.5_f32 * 0.5 + 0.5 * 0.5).sqrt();
        assert!((c.radius - expected_r).abs() < 1e-4);
    }

    #[test]
    fn test_equilateral_triangle() {
        let h = (3.0_f32).sqrt() / 2.0;
        let c = smallest_enclosing_circle(&[
            Point(0.0, 0.0),
            Point(1.0, 0.0),
            Point(0.5, h),
        ]);
        // All three points should be on the circle
        let d1 = c.center.distance_to(&Point(0.0, 0.0));
        let d2 = c.center.distance_to(&Point(1.0, 0.0));
        let d3 = c.center.distance_to(&Point(0.5, h));
        assert!((d1 - c.radius).abs() < 1e-4);
        assert!((d2 - c.radius).abs() < 1e-4);
        assert!((d3 - c.radius).abs() < 1e-4);
    }

    #[test]
    fn test_all_points_contained() {
        let points = vec![
            Point(0.0, 0.0),
            Point(3.0, 1.0),
            Point(1.0, 4.0),
            Point(-1.0, 2.0),
            Point(2.0, 2.0),
        ];
        let c = smallest_enclosing_circle(&points);
        for p in &points {
            let d = c.center.distance_to(p);
            assert!(d <= c.radius + 1e-4, "Point {:?} outside circle (d={}, r={})", p, d, c.radius);
        }
    }
}
