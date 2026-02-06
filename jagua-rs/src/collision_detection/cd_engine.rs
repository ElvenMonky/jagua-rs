use crate::collision_detection::hazards::collector::HazardCollector;
use crate::collision_detection::hazards::filter::HazardFilter;
use crate::collision_detection::hazards::{HazKey, Hazard, HazardEntity};
use crate::entities::PItemKey;
use crate::geometry::Transformation;
use crate::geometry::fail_fast::{SPSurrogate, SPSurrogateConfig};
use crate::geometry::geo_enums::{GeoPosition, GeoRelation};
use crate::geometry::geo_traits::{CollidesWith, Transformable};
use crate::geometry::primitives::{Circle, Edge, Rect, SPolygon};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use slotmap::SlotMap;
use std::collections::HashSet;

// ─── Axis sweep ─────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum AxisMark {
    Min = 0,
    Max = 1,
}

#[derive(Clone, Copy, Debug)]
struct AxisEvent {
    value: f32,
    mark: AxisMark,
    hkey: HazKey,
}

impl AxisEvent {
    #[inline]
    fn sort_key(&self) -> (ordered_float::OrderedFloat<f32>, AxisMark) {
        (ordered_float::OrderedFloat(self.value), self.mark)
    }
}

#[derive(Clone, Debug, Default)]
struct AxisList {
    events: Vec<AxisEvent>,
}

impl AxisList {
    fn new() -> Self {
        Self::default()
    }

    fn insert(&mut self, hkey: HazKey, min_val: f32, max_val: f32) {
        let min_event = AxisEvent { value: min_val, mark: AxisMark::Min, hkey };
        let max_event = AxisEvent { value: max_val, mark: AxisMark::Max, hkey };
        let min_pos = self.insertion_point(&min_event);
        self.events.insert(min_pos, min_event);
        let max_pos = self.insertion_point(&max_event);
        self.events.insert(max_pos, max_event);
    }

    fn remove(&mut self, hkey: HazKey) {
        self.events.retain(|e| e.hkey != hkey);
    }

    #[inline]
    fn insertion_point(&self, event: &AxisEvent) -> usize {
        let key = event.sort_key();
        self.events.partition_point(|e| e.sort_key() < key)
    }

    /// Single-pass sweep producing three candidate sets for this axis:
    ///
    /// - `all_overlap`: all hazards whose interval overlaps [e_min, e_max]
    /// - `surrounds_entity`: hazards whose interval fully surrounds [e_min, e_max]
    /// - `inside_entity`: hazards whose interval is fully inside [e_min, e_max]
    fn compute_candidates(
        &self,
        e_min: f32,
        e_max: f32,
    ) -> (HashSet<HazKey>, HashSet<HazKey>, HashSet<HazKey>) {
        let mut all_overlap = HashSet::new();
        let mut surrounds_entity = HashSet::new();
        let mut inside_entity = HashSet::new();

        for event in &self.events {
            if event.value > e_max {
                break;
            }
            match (event.mark, event.value <= e_min) {
                (AxisMark::Min, true) => {
                    all_overlap.insert(event.hkey);
                    surrounds_entity.insert(event.hkey);
                }
                (AxisMark::Min, false) => {
                    all_overlap.insert(event.hkey);
                }
                (AxisMark::Max, true) => {
                    all_overlap.remove(&event.hkey);
                    surrounds_entity.remove(&event.hkey);
                }
                (AxisMark::Max, false) => {
                    if surrounds_entity.contains(&event.hkey) {
                        surrounds_entity.remove(&event.hkey);
                    } else {
                        inside_entity.insert(event.hkey);
                    }
                }
            }
        }

        (all_overlap, surrounds_entity, inside_entity)
    }
}

// ─── Generic hazard collision: surrogate-first ──────────────────

/// Check if entity collides with hazard shape via edge/surrogate intersection.
/// Check order: hazard bbox → hazard surrogate → hazard edges.
#[inline]
fn entity_collides_with_hazard<T>(entity: &T, haz_shape: &SPolygon) -> bool
where
    T: CollidesWith<Circle> + CollidesWith<Edge> + CollidesWith<Rect>,
{
    if !entity.collides_with(&haz_shape.bbox) {
        return false;
    }
    if let Some(surr) = &haz_shape.surrogate {
        for pole in surr.ff_poles() {
            if entity.collides_with(pole) {
                return true;
            }
        }
        for pier in surr.ff_piers() {
            if entity.collides_with(pier) {
                return true;
            }
        }
    }
    for haz_edge in haz_shape.edge_iter() {
        if entity.collides_with(&haz_edge) {
            return true;
        }
    }
    false
}

// ─── Containment checks ────────────────────────────────────────

/// Check if entity is fully contained inside hazard.
/// Only valid when hazard bbox surrounds entity bbox.
#[inline]
fn entity_contained_in_hazard(
    entity_poi_center: &crate::geometry::primitives::Point,
    haz_shape: &SPolygon,
    haz_entity: HazardEntity,
) -> bool {
    let contained = haz_shape.collides_with(entity_poi_center);
    match (haz_entity.scope(), contained) {
        (GeoPosition::Interior, true) | (GeoPosition::Exterior, false) => true,
        _ => false,
    }
}

/// Check if hazard is fully contained inside entity.
/// Only valid when entity bbox surrounds hazard bbox.
#[inline]
fn hazard_contained_in_entity(
    entity: &SPolygon,
    haz_shape: &SPolygon,
    haz_entity: HazardEntity,
) -> bool {
    let contained = entity.collides_with(&haz_shape.poi.center);
    match (haz_entity.scope(), contained) {
        (GeoPosition::Interior, true) | (GeoPosition::Exterior, false) => true,
        _ => false,
    }
}

// ─── CDEngine ──────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct CDEngine {
    x_axis: AxisList,
    y_axis: AxisList,
    pub hazards_map: SlotMap<HazKey, Hazard>,
    pub config: CDEConfig,
    bbox: Rect,
    hkey_exterior: HazKey,
}

impl CDEngine {
    pub fn new(bbox: Rect, static_hazards: Vec<Hazard>, config: CDEConfig) -> Self {
        let mut engine = CDEngine {
            x_axis: AxisList::new(),
            y_axis: AxisList::new(),
            hazards_map: SlotMap::with_key(),
            config,
            bbox,
            hkey_exterior: HazKey::default(),
        };

        for haz in static_hazards {
            engine.register_hazard(haz);
        }

        engine.hkey_exterior = engine
            .hazards_map
            .iter()
            .find(|(_, h)| matches!(h.entity, HazardEntity::Exterior))
            .map(|(hkey, _)| hkey)
            .expect("No exterior hazard registered");

        engine
    }

    pub fn register_hazard(&mut self, hazard: Hazard) {
        debug_assert!(
            !self.hazards_map.values().any(|h| h.entity == hazard.entity),
            "Hazard with an identical entity already registered"
        );
        let haz_bbox = hazard.shape.bbox;
        let hkey = self.hazards_map.insert(hazard);
        self.x_axis.insert(hkey, haz_bbox.x_min, haz_bbox.x_max);
        self.y_axis.insert(hkey, haz_bbox.y_min, haz_bbox.y_max);
    }

    pub fn deregister_hazard_by_entity(&mut self, hazard_entity: HazardEntity) -> Hazard {
        let hkey = self
            .hazards_map
            .iter()
            .find(|(_, h)| h.entity == hazard_entity)
            .map(|(hkey, _)| hkey)
            .expect("Cannot deregister hazard that is not registered");
        self.deregister_hazard_by_key(hkey)
    }

    pub fn deregister_hazard_by_key(&mut self, hkey: HazKey) -> Hazard {
        let hazard = self
            .hazards_map
            .remove(hkey)
            .expect("Cannot deregister hazard that is not registered");
        self.x_axis.remove(hkey);
        self.y_axis.remove(hkey);
        hazard
    }

    pub fn save(&self) -> CDESnapshot {
        let dynamic_hazards = self
            .hazards_map
            .values()
            .filter(|h| h.dynamic)
            .cloned()
            .collect_vec();
        CDESnapshot { dynamic_hazards }
    }

    pub fn restore(&mut self, snapshot: &CDESnapshot) {
        let mut hazards_to_remove = self
            .hazards_map
            .iter()
            .filter(|(_, h)| h.dynamic)
            .map(|(hkey, h)| (hkey, h.entity))
            .collect_vec();
        let mut hazards_to_add = vec![];

        for hazard in snapshot.dynamic_hazards.iter() {
            let present = hazards_to_remove
                .iter()
                .position(|(_, h)| h == &hazard.entity);
            if let Some(idx) = present {
                hazards_to_remove.swap_remove(idx);
            } else {
                hazards_to_add.push(hazard.clone());
            }
        }

        for (hkey, _) in hazards_to_remove {
            self.deregister_hazard_by_key(hkey);
        }
        for hazard in hazards_to_add {
            self.register_hazard(hazard);
        }
    }

    pub fn hazards(&self) -> impl Iterator<Item = &Hazard> {
        self.hazards_map.values()
    }

    // ─── Candidate computation ──────────────────────────────────

    /// Compute collision candidates by intersecting per-axis results.
    /// Returns (edge_intersection, entity_inside_hazard, hazard_inside_entity).
    fn compute_candidates(
        &self,
        entity_bbox: Rect,
    ) -> (HashSet<HazKey>, HashSet<HazKey>, HashSet<HazKey>) {
        let (x_all, x_surr, x_inside) =
            self.x_axis.compute_candidates(entity_bbox.x_min, entity_bbox.x_max);
        let (y_all, y_surr, y_inside) =
            self.y_axis.compute_candidates(entity_bbox.y_min, entity_bbox.y_max);

        let edge_intersection: HashSet<HazKey> =
            x_all.intersection(&y_all).copied().collect();
        let entity_inside_hazard: HashSet<HazKey> =
            x_surr.intersection(&y_surr).copied().collect();
        let hazard_inside_entity: HashSet<HazKey> =
            x_inside.intersection(&y_inside).copied().collect();

        (edge_intersection, entity_inside_hazard, hazard_inside_entity)
    }

    // ─── Surrogate bbox ─────────────────────────────────────────

    fn compute_surrogate_bbox(
        base_surrogate: &SPSurrogate,
        transform: &Transformation,
    ) -> Rect {
        let mut x_min = f32::MAX;
        let mut y_min = f32::MAX;
        let mut x_max = f32::MIN;
        let mut y_max = f32::MIN;

        for pole in base_surrogate.ff_poles() {
            let b = pole.transform_clone(transform).bbox();
            x_min = x_min.min(b.x_min);
            y_min = y_min.min(b.y_min);
            x_max = x_max.max(b.x_max);
            y_max = y_max.max(b.y_max);
        }
        for pier in base_surrogate.ff_piers() {
            let b = pier.transform_clone(transform).bbox();
            x_min = x_min.min(b.x_min);
            y_min = y_min.min(b.y_min);
            x_max = x_max.max(b.x_max);
            y_max = y_max.max(b.y_max);
        }

        Rect { x_min, y_min, x_max, y_max }
    }

    pub fn detect_poly_collision(&self, shape: &SPolygon, filter: &impl HazardFilter) -> bool {
        if self.bbox.relation_to(shape.bbox) != GeoRelation::Surrounding {
            return true;
        }

        let (edge_intersection, entity_inside_hazard, hazard_inside_entity) =
            self.compute_candidates(shape.bbox);

        for &hkey in &edge_intersection {
            if filter.is_irrelevant(hkey) {
                continue;
            }
            let hazard = &self.hazards_map[hkey];

            if entity_collides_with_hazard(shape, &hazard.shape) {
                return true;
            }
            if entity_inside_hazard.contains(&hkey)
                && entity_contained_in_hazard(&shape.poi.center, &hazard.shape, hazard.entity)
            {
                return true;
            }
            if hazard_inside_entity.contains(&hkey)
                && hazard_contained_in_entity(shape, &hazard.shape, hazard.entity)
            {
                return true;
            }
        }

        false
    }

    // ─── Phase 1: surrogate as entity ───────────────────────────
    //
    // No containment checks — surrogate doesn't represent full shape.
    // All candidates get intersection check only.

    pub fn detect_surrogate_collision(
        &self,
        base_surrogate: &SPSurrogate,
        transform: &Transformation,
        filter: &impl HazardFilter,
    ) -> bool {
        let surrogate_bbox = Self::compute_surrogate_bbox(base_surrogate, transform);
        let (edge_intersection, _, _) = self.compute_candidates(surrogate_bbox);

        for &hkey in &edge_intersection {
            if filter.is_irrelevant(hkey) {
                continue;
            }
            let haz_shape = &self.hazards_map[hkey].shape;
            for pole in base_surrogate.ff_poles() {
                let t_pole = pole.transform_clone(transform);
                if entity_collides_with_hazard(&t_pole, haz_shape) {
                    return true;
                }
            }
            for pier in base_surrogate.ff_piers() {
                let t_pier = pier.transform_clone(transform);
                if entity_collides_with_hazard(&t_pier, haz_shape) {
                    return true;
                }
            }
        }

        false
    }

    // ─── Phase 2: full polygon as entity ────────────────────────
    //
    // Entity's own surrogate is NOT re-checked (phase 1 already did that).
    // Hazard surrogate IS used via entity_collides_with_hazard.
    //
    // Per-category checks:
    //   edge_intersection    → intersection only
    //   entity_inside_hazard → intersection + entity_contained_in_hazard
    //   hazard_inside_entity → intersection + hazard_contained_in_entity

    pub fn collect_poly_collisions(&self, shape: &SPolygon, collector: &mut impl HazardCollector) {
        if self.bbox.relation_to(shape.bbox) != GeoRelation::Surrounding {
            collector.insert(self.hkey_exterior, HazardEntity::Exterior);
        }

        let (edge_intersection, entity_inside_hazard, hazard_inside_entity) =
            self.compute_candidates(shape.bbox);

        for &hkey in &edge_intersection {
            if collector.contains_key(hkey) {
                continue;
            }
            let hazard = &self.hazards_map[hkey];

            if entity_collides_with_hazard(shape, &hazard.shape)
                || (entity_inside_hazard.contains(&hkey)
                    && entity_contained_in_hazard(&shape.poi.center, &hazard.shape, hazard.entity))
                || (hazard_inside_entity.contains(&hkey)
                    && hazard_contained_in_entity(shape, &hazard.shape, hazard.entity))
            {
                collector.insert(hkey, hazard.entity);
            }
        }
    }

    pub fn collect_surrogate_collisions(
        &self,
        base_surrogate: &SPSurrogate,
        transform: &Transformation,
        collector: &mut impl HazardCollector,
    ) {
        let surrogate_bbox = Self::compute_surrogate_bbox(base_surrogate, transform);
        let (edge_intersection, _, _) = self.compute_candidates(surrogate_bbox);

        for &hkey in &edge_intersection {
            if collector.contains_key(hkey) {
                continue;
            }
            let haz_shape = &self.hazards_map[hkey].shape;
            let mut collides = false;
            for pole in base_surrogate.ff_poles() {
                if collides { break; }
                let t_pole = pole.transform_clone(transform);
                collides = entity_collides_with_hazard(&t_pole, haz_shape);
            }
            for pier in base_surrogate.ff_piers() {
                if collides { break; }
                let t_pier = pier.transform_clone(transform);
                collides = entity_collides_with_hazard(&t_pier, haz_shape);
            }
            if collides {
                collector.insert(hkey, self.hazards_map[hkey].entity);
            }
        }
    }

        pub fn bbox(&self) -> Rect {
        self.bbox
    }

    pub fn haz_key_from_pi_key(&self, pik: PItemKey) -> Option<HazKey> {
        self.hazards_map
            .iter()
            .find(|(_, hazard)| match hazard.entity {
                HazardEntity::PlacedItem { pk, .. } => pik == pk,
                _ => false,
            })
            .map(|(key, _)| key)
    }
}

///Configuration of the [`CDEngine`]
#[derive(Serialize, Deserialize, Clone, Copy, Debug, PartialEq)]
pub struct CDEConfig {
    ///Maximum depth of the quadtree
    pub quadtree_depth: u8,
    /// Stop traversing the quadtree and perform collision collection immediately when the total number of edges in a node falls below this number
    pub cd_threshold: u8,
    ///Configuration of the surrogate generation for items
    pub item_surrogate_config: SPSurrogateConfig,
}

/// Snapshot of the state of [`CDEngine`]. Can be used to restore to a previous state.
#[derive(Clone, Debug)]
pub struct CDESnapshot {
    pub dynamic_hazards: Vec<Hazard>,
}
