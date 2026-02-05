use crate::collision_detection::hazards::HazKey;
use crate::collision_detection::hazards::filter::HazardFilter;
use crate::collision_detection::quadtree::QTHazPresence;
use crate::collision_detection::quadtree::QTHazard;
use std::cmp::Ordering;
use std::ops::Not;

/// Vector of `QTHazard`s, which always remains sorted by presence.
/// <br>
/// This is a performance optimization to be able to quickly return the "strongest" hazard
/// Strongest meaning the highest [`QTHazPresence`] (`Entire` > `Partial(high)` > `Partial(low)` > `None`)
#[derive(Clone, Debug, Default)]
pub struct QTHazardVec {
    hazards: Vec<QTHazard>,
    /// Number of edges from active hazards in the vector
    n_active_edges: usize,
}

impl QTHazardVec {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, haz: QTHazard) {
        debug_assert!(!matches!(haz.presence, QTHazPresence::None));
        debug_assert!(
            self.hazards
                .iter()
                .filter(|other| other.entity == haz.entity || other.hkey == haz.hkey)
                .count()
                == 0,
            "More than one hazard from same entity or key in the vector! (This should never happen!)"
        );
        match self
            .hazards
            .binary_search_by(|probe| order_by_descending_strength(probe, &haz))
        {
            Ok(pos) | Err(pos) => {
                self.n_active_edges += haz.n_edges();
                self.hazards.insert(pos, haz);
            }
        }
    }

    pub fn remove(&mut self, hkey: HazKey) -> Option<QTHazard> {
        let pos = self.hazards.iter().position(|ch| ch.hkey == hkey);
        match pos {
            Some(pos) => {
                let haz = self.hazards.remove(pos);
                self.n_active_edges -= haz.n_edges();
                Some(haz)
            }
            None => None,
        }
    }

    #[inline(always)]
    /// Returns the strongest hazard (if any) (`Entire` > `Partial(high)` > `Partial(low)` > `None`)
    /// Ignores any hazards that are deemed irrelevant by the filter.
    pub fn strongest(&self, filter: &impl HazardFilter) -> Option<&QTHazard> {
        debug_assert!(assert_caches_correct(self));
        self.iter().find(|hz| !filter.is_irrelevant(hz.hkey))
    }

    pub fn first(&self) -> Option<&QTHazard> {
        self.hazards.first()
    }

    pub fn is_empty(&self) -> bool {
        self.hazards.is_empty()
    }

    pub fn len(&self) -> usize {
        self.hazards.len()
    }
    pub fn no_partial_hazards(&self) -> bool {
        self.hazards
            .iter()
            .any(|hz| matches!(hz.presence, QTHazPresence::Partial(_)))
            .not()
    }

    pub fn iter(&self) -> impl Iterator<Item = &QTHazard> {
        self.hazards.iter()
    }

    pub fn n_active_edges(&self) -> usize {
        debug_assert!(assert_caches_correct(self));
        self.n_active_edges
    }
}

fn order_by_descending_strength(qth1: &QTHazard, qth2: &QTHazard) -> Ordering {
    // Sort by presence: Entire (1.0 + full area) > Partial(high presence area) > Partial(low presence area) > None (-1.0)
    let qth_presence_sortkey = |qth: &QTHazard| -> ordered_float::OrderedFloat<f32> {
        match &qth.presence {
            QTHazPresence::None => ordered_float::OrderedFloat(-1.0),
            QTHazPresence::Partial(p) => ordered_float::OrderedFloat(p.presence_area),
            QTHazPresence::Entire => ordered_float::OrderedFloat(1.0 + qth.qt_bbox.area()),
        }
    };

    qth_presence_sortkey(qth1)
        .cmp(&qth_presence_sortkey(qth2))
        .reverse()
}

fn assert_caches_correct(qthazard_vec: &QTHazardVec) -> bool {
    assert!(
        qthazard_vec
            .hazards
            .windows(2)
            .all(|w| order_by_descending_strength(&w[0], &w[1]) != Ordering::Greater),
        "Hazards are not sorted correctly!"
    );
    assert_eq!(
        qthazard_vec
            .hazards
            .iter()
            .map(|hz| hz.n_edges())
            .sum::<usize>(),
        qthazard_vec.n_active_edges,
        "Active edges count is not correct!"
    );
    true
}
