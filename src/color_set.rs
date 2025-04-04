use std::collections::HashSet;
use std::fmt;

const BITS_COUNT: i32 = 22;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct ColorSet(HashSet<u64>);

impl ColorSet {
    pub fn from(value: u64, prev_count: usize, curr_count: usize) -> ColorSet {
        let mut result = ColorSet(HashSet::new());
        result.insert(value, prev_count, curr_count);
        result
    }

    pub fn equal_set(&self, other: &ColorSet) -> bool {
        self.0.eq(&other.0)
    }

    pub fn intersection(&self, other: &ColorSet, _intersect: u8) -> ColorSet {
        let sub0 = self.sub_intersection(other, 0);
        let sub1 = self.sub_intersection(other, 1);
        let result = sub0.0.union(&sub1.0).copied().collect();
        ColorSet(result)
    }

    fn sub_intersection(&self, other: &ColorSet, intersect: u8) -> ColorSet {
        let keys: HashSet<u64> = other
            .0
            .iter()
            .map(|p| match intersect {
                0 => !((((1 << BITS_COUNT) - 1) as u64) << BITS_COUNT) & p,
                1 => !(((1 << BITS_COUNT) - 1) as u64) & p,
                _ => unreachable!("Should never be reached"),
            })
            .collect();
        let result = self
            .0
            .iter()
            .filter(|p| match intersect {
                0 => keys.contains(&(!((((1 << BITS_COUNT) - 1) as u64) << BITS_COUNT) & **p)),
                1 => keys.contains(&(!(((1 << BITS_COUNT) - 1) as u64) & **p)),
                _ => unreachable!("Should never be reached"),
            })
            .copied()
            .collect();
        ColorSet(result)
    }

    pub fn true_intersection(&self, other: &ColorSet) -> ColorSet {
        ColorSet(self.0.intersection(&other.0).copied().collect())
    }

    pub fn difference(&self, other: &ColorSet) -> ColorSet {
        ColorSet(self.0.difference(&other.0).copied().collect::<HashSet<_>>())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn insert(&mut self, value: u64, prev_count: usize, curr_count: usize) {
        if prev_count >= (1 << BITS_COUNT) {
            log::error!("prev_count is {}", prev_count);
        }
        assert!(prev_count < (1 << BITS_COUNT));
        if curr_count >= (1 << BITS_COUNT) {
            log::error!("curr_count is {}", curr_count);
        }
        assert!(curr_count < (1 << BITS_COUNT));
        let value =
            (value << (2 * BITS_COUNT)) | ((prev_count as u64) << BITS_COUNT) | curr_count as u64;
        self.0.insert(value);
    }

    pub fn contains(&self, value: u64, prev_count: usize, curr_count: usize) -> bool {
        if prev_count >= (1 << BITS_COUNT) {
            log::error!("prev_count is {}", prev_count);
        }
        assert!(prev_count < (1 << BITS_COUNT));
        if curr_count >= (1 << BITS_COUNT) {
            log::error!("curr_count is {}", curr_count);
        }
        assert!(curr_count < (1 << BITS_COUNT));
        let value =
            (value << (2 * BITS_COUNT)) | ((prev_count as u64) << BITS_COUNT) | curr_count as u64;
        self.0.contains(&value)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn new() -> Self {
        Self(HashSet::new())
    }
}

impl PartialOrd for ColorSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ColorSet {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.len().cmp(&other.0.len())
    }
}

impl fmt::Display for ColorSet {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{ ")?;
        for el in self.0.iter() {
            write!(
                f,
                "{}|{}|{}, ",
                el >> (2 * BITS_COUNT),
                (el >> BITS_COUNT) & ((1 << BITS_COUNT) - 1),
                (el) & ((1 << BITS_COUNT) - 1)
            )?;
        }
        write!(f, "}} ")
    }
}
