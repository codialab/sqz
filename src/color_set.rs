use std::collections::HashSet;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct ColorSet(HashSet<u64>);

impl ColorSet {
    pub fn from(value: u64, prev_count: usize, curr_count: usize) -> ColorSet {
        let mut result = ColorSet(HashSet::new());
        result.insert(value, prev_count, curr_count);
        result
    }

    pub fn intersection(&self, other: &ColorSet, intersect: u8) -> ColorSet {
        let keys: HashSet<u64> = other
            .0
            .iter()
            .map(|p| match intersect {
                0 => !(1048575u64 << 20) & p,
                1 => !1048575u64 & p,
                _ => unreachable!("Should never be reached"),
            })
            .collect();
        let result = self
            .0
            .iter()
            .filter(|p| match intersect {
                0 => keys.contains(&(!(1048575u64 << 20) & **p)),
                1 => keys.contains(&(!1048575u64 & **p)),
                _ => unreachable!("Should never be reached"),
            })
            .copied()
            .collect();
        ColorSet(result)
    }

    pub fn difference(&self, other: &ColorSet) -> ColorSet {
        ColorSet(self.0.difference(&other.0).copied().collect::<HashSet<_>>())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn insert(&mut self, value: u64, prev_count: usize, curr_count: usize) {
        assert!(prev_count < 1048576);
        assert!(curr_count < 1048576);
        let value = (value << 40) | ((prev_count as u64) << 20) | curr_count as u64;
        self.0.insert(value);
    }

    pub fn len(&self) -> usize {
        self.0.len()
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
