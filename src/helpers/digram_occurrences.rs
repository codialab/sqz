use std::{
    collections::{hash_map::Iter, HashSet},
    hash::{BuildHasherDefault, DefaultHasher},
};

use anyhow::{anyhow, Result};
use deepsize::DeepSizeOf;

use crate::helpers::{
    utils::LocalizedDigram, AddressNumber, CanonicalDigram, DeterministicHashMap,
    DeterministicHashSet, Freq, Neighbors, NodeId, Occurrence, Side,
};

#[derive(Clone, Debug, DeepSizeOf)]
pub struct DigramOccurrences {
    inner: DeterministicHashMap<CanonicalDigram, DeterministicHashSet<Occurrence>>,
    neighbors: Neighbors,
    freq: Freq,
}

impl DigramOccurrences {
    pub fn new() -> Self {
        Self {
            inner: DeterministicHashMap::default(),
            neighbors: Neighbors::new(),
            freq: Freq::new(),
        }
    }

    pub fn print_occurrences(&self) {
        for (digram, occurrences) in &self.inner {
            println!("{:?}: {:?}", digram, occurrences);
        }
    }

    pub fn total_len(&self) -> usize {
        self.into_iter().map(|(_, x)| x.len()).sum()
    }

    #[allow(dead_code)]
    pub fn contains(&self, k: &CanonicalDigram) -> bool {
        self.inner.contains_key(k)
    }

    #[allow(dead_code)]
    pub fn contains_occurrence(&self, k: &CanonicalDigram, o: &Occurrence) -> bool {
        self.contains(k) && self.inner[k].contains(o)
    }

    pub fn shrink(&mut self) {
        self.inner.shrink_to_fit();
        self.neighbors.shrink();
        self.freq.shrink();
    }

    pub fn log_sizes(&self) {
        log::info!(
            "DigramOccurrences with {} occurrences (inner: {}MB, neighbors: {}MB, freq: {}MB)",
            self.total_len(),
            self.inner.deep_size_of() / 1_000_000,
            self.neighbors.deep_size_of() / 1_000_000,
            self.freq.deep_size_of() / 1_000_000
        );
    }

    pub fn from(haplotypes: Vec<Vec<LocalizedDigram>>) -> Self {
        let mut d = Self::new();

        for (i, haplotype) in haplotypes.into_iter().enumerate() {
            for local_digram in haplotype {
                let (digram, address) = local_digram.split_to_canonical();
                d.insert(&digram, Occurrence::new(i, address));
                d.neighbors.insert(
                    Side::Right,
                    (digram.get_u(), i, address.get_first()),
                    (digram.get_v(), address.get_second()),
                );
                d.neighbors.insert(
                    Side::Left,
                    (digram.get_v(), i, address.get_second()),
                    (digram.get_u(), address.get_first()),
                );
            }
        }

        d.freq = d.get_freq();
        d.shrink();
        d.log_sizes();
        d
    }

    pub fn remove_digram(
        &mut self,
        index: &CanonicalDigram,
    ) -> HashSet<Occurrence, BuildHasherDefault<DefaultHasher>> {
        let (_, occurrences) = self
            .inner
            .remove_entry(index)
            .expect("DigramOccurrences needs to contain a digram to remove it");
        let frequency = occurrences.len();
        self.freq.remove(index, frequency);
        self.neighbors
            .remove_occurrences(Side::Left, index, &occurrences);
        self.neighbors
            .remove_occurrences(Side::Right, index, &occurrences);
        occurrences
    }

    #[cfg(test)]
    pub fn are_neighbors_equal(&self) -> bool {
        self.neighbors.inner.iter().all(|(_, (v1, v2))| v1 == v2)
    }

    pub fn delete_occurrence(&mut self, index: &CanonicalDigram, value: &Occurrence) -> Result<()> {
        let old_occurrence = self.inner[index].len();
        let mut has_deleted = false;
        if index.is_symmetric() {
            self.neighbors
                .remove_occurrence(Side::Left, index, &value.flip());
            self.neighbors
                .remove_occurrence(Side::Right, index, &value.flip());
            if self
                .inner
                .get_mut(index)
                .expect(
                    "DigramOccurrences needs to contain digram to remove one of its occurrences",
                )
                .remove(&value.flip())
            {
                has_deleted = true;
            }
        }
        self.neighbors.remove_occurrence(Side::Left, index, value);
        self.neighbors.remove_occurrence(Side::Right, index, value);
        self.freq
            .change_frequency(index, old_occurrence, old_occurrence - 1);
        if self
            .inner
            .get_mut(index)
            .expect("DigramOccurrences needs to contain digram to remove one of its occurrences")
            .remove(value)
        {
            has_deleted = true;
        }
        if !has_deleted {
            return Err(anyhow!(
                "Cannot delete occurrence in DigramOccurrences.inner"
            ));
        }
        if self.inner[index].is_empty() {
            return match self.inner.remove(index) {
                Some(_) => Ok(()),
                None => Err(anyhow!(
                    "Cannot delete occurrence set in DigramOccurrences.inner"
                )),
            };
        }
        Ok(())
    }

    // TODO change freq to not contain digrams of frequency 0

    pub fn add_digram(
        &mut self,
        index: &CanonicalDigram,
        values: HashSet<Occurrence, BuildHasherDefault<DefaultHasher>>,
    ) {
        if values.is_empty() {
            return;
        }
        self.freq.insert(index.clone(), values.len());
        self.neighbors
            .insert_values(Side::Left, index, values.iter().cloned().collect());
        self.neighbors
            .insert_values(Side::Right, index, values.iter().cloned().collect());
        self.inner.insert(index.clone(), values);
    }

    pub fn add_occurrence(&mut self, index: &CanonicalDigram, value: Occurrence) {
        self.neighbors
            .insert_occurrence(Side::Left, index, value.clone());
        self.neighbors
            .insert_occurrence(Side::Right, index, value.clone());
        if self.inner.contains_key(index) {
            let old_occurrence = self.inner[index].len();
            self.freq
                .change_frequency(index, old_occurrence, old_occurrence + 1);
        } else {
            self.freq.insert(index.clone(), 1);
        }
        self.inner.entry(index.clone()).or_default().insert(value);
    }

    pub fn get_left_neighbor(
        &self,
        index: &CanonicalDigram,
        value: &Occurrence,
    ) -> Option<(NodeId, AddressNumber)> {
        let node_id = index.get_u();
        let haplotype_id = value.0;
        let address_number = value.1 .0;
        self.neighbors
            .get_value(Side::Left, node_id, haplotype_id as usize, address_number)
    }

    pub fn get_right_neighbor(
        &self,
        index: &CanonicalDigram,
        value: &Occurrence,
    ) -> Option<(NodeId, AddressNumber)> {
        let node_id = index.get_v();
        let haplotype_id = value.0;
        let address_number = value.1 .1;
        self.neighbors
            .get_value(Side::Right, node_id, haplotype_id as usize, address_number)
    }

    pub fn get_most_frequent(&self, minimum_freq: usize) -> Option<CanonicalDigram> {
        let most_freq = self.freq.get_most_freq()?;
        if self.inner[&most_freq].len() < minimum_freq {
            None
        } else {
            Some(most_freq)
        }
    }

    fn insert(&mut self, index: &CanonicalDigram, value: Occurrence) {
        self.inner.entry(index.clone()).or_default().insert(value);
    }

    fn get_freq(&self) -> Freq {
        let mut freq = Freq::new();
        for (digram, occurrences) in &self.inner {
            freq.insert(digram.clone(), occurrences.len());
        }
        freq
    }
}

pub struct DigramOccurrencesRefIter<'a> {
    iter: Iter<'a, CanonicalDigram, DeterministicHashSet<Occurrence>>,
}

impl<'a> Iterator for DigramOccurrencesRefIter<'a> {
    // HashMap iterators yield a tuple: (&Key, &Value)
    type Item = (&'a CanonicalDigram, &'a DeterministicHashSet<Occurrence>);

    fn next(&mut self) -> Option<Self::Item> {
        // We simply delegate the work to the inner iterator!
        self.iter.next()
    }
}

impl<'a> IntoIterator for &'a DigramOccurrences {
    type Item = (&'a CanonicalDigram, &'a DeterministicHashSet<Occurrence>);
    type IntoIter = DigramOccurrencesRefIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        DigramOccurrencesRefIter {
            iter: self.inner.iter(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::helpers::utils::{Address, Digram, LocalizedDigram, Orientation, UndirectedNodeId};

    use super::*;

    fn get_canonical_digram() -> CanonicalDigram {
        get_digram().into()
    }

    fn get_digram() -> Digram {
        Digram(
            NodeId(UndirectedNodeId::new(1), Orientation::Forward),
            NodeId(UndirectedNodeId::new(2), Orientation::Forward),
        )
    }

    fn get_doc() -> DigramOccurrences {
        let doc = DigramOccurrences::from(vec![vec![LocalizedDigram(
            get_digram(),
            Address(AddressNumber(0), AddressNumber(1)),
        )]]);
        doc
    }

    fn get_digram_from(a: u32, b: u32) -> Digram {
        Digram(
            NodeId(UndirectedNodeId::new(a), Orientation::Forward),
            NodeId(UndirectedNodeId::new(b), Orientation::Forward),
        )
    }

    #[test]
    fn test_insert() {
        let mut doc = DigramOccurrences::new();
        assert_eq!(doc.inner.len(), 0);
        doc.insert(
            &get_canonical_digram(),
            Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        );
        assert_eq!(doc.inner.len(), 1);
        assert_eq!(doc.neighbors.left_len(), 0);
        assert_eq!(doc.neighbors.right_len(), 0);
        assert_eq!(doc.freq.inner.len(), 0);
    }

    #[test]
    fn test_contains() {
        let mut doc = DigramOccurrences::new();
        assert!(!doc.contains(&get_canonical_digram()));
        doc.insert(
            &get_canonical_digram(),
            Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        );
        assert!(doc.contains(&get_canonical_digram()));
    }

    #[test]
    fn test_from() {
        let doc = DigramOccurrences::from(vec![vec![LocalizedDigram(
            get_digram(),
            Address(AddressNumber(0), AddressNumber(1)),
        )]]);
        assert_eq!(doc.inner.len(), 1);
        assert_eq!(doc.neighbors.left_len(), 1);
        assert_eq!(doc.neighbors.right_len(), 1);
        assert_eq!(doc.freq.inner.len(), 1);
    }

    #[test]
    fn test_remove_digram() {
        let mut doc = get_doc();
        doc.remove_digram(&get_canonical_digram());
        assert_eq!(doc.inner.len(), 0);
        assert_eq!(doc.neighbors.left_len(), 0);
        assert_eq!(doc.neighbors.right_len(), 0);
        assert_eq!(doc.freq.inner.len(), 0);
    }

    #[test]
    fn test_delete_occurrence() {
        let mut doc = get_doc();
        doc.delete_occurrence(
            &get_canonical_digram(),
            &Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        )
        .expect("Can delete occurrence");
        assert_eq!(doc.inner.len(), 0);
        assert_eq!(doc.neighbors.left_len(), 0);
        assert_eq!(doc.neighbors.right_len(), 0);
        assert_eq!(doc.freq.inner.len(), 0);
    }

    #[test]
    fn test_add_digram() {
        let mut doc = DigramOccurrences::new();
        let mut occurrences: DeterministicHashSet<Occurrence> = DeterministicHashSet::default();
        occurrences.insert(Occurrence(0, Address(AddressNumber(0), AddressNumber(1))));
        doc.add_digram(&get_canonical_digram(), occurrences);
        assert_eq!(doc.inner.len(), 1);
        assert_eq!(doc.neighbors.left_len(), 1);
        assert_eq!(doc.neighbors.right_len(), 1);
        assert_eq!(doc.freq.inner.len(), 1);
    }

    #[test]
    fn test_add_occurrence() {
        let mut doc = DigramOccurrences::new();
        doc.add_occurrence(
            &get_canonical_digram(),
            Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        );
        assert_eq!(doc.inner.len(), 1);
        assert_eq!(doc.neighbors.left_len(), 1);
        assert_eq!(doc.neighbors.right_len(), 1);
        assert_eq!(doc.freq.inner.len(), 1);
    }

    #[test]
    fn test_get_left_neighbor() {
        let mut doc = get_doc();
        doc.add_occurrence(
            &get_digram_from(0, 1).into(),
            Occurrence(0, Address(AddressNumber(3), AddressNumber(0))),
        );
        let ln = doc.get_left_neighbor(
            &get_canonical_digram(),
            &Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        );
        assert!(ln.is_some());
        let (node_id, address_number) = ln.unwrap();
        assert_eq!(
            node_id,
            NodeId(UndirectedNodeId::new(0), Orientation::Forward)
        );
        assert_eq!(address_number, AddressNumber(3));
    }

    #[test]
    fn test_get_right_neighbor() {
        let mut doc = get_doc();
        doc.add_occurrence(
            &get_digram_from(2, 3).into(),
            Occurrence(0, Address(AddressNumber(1), AddressNumber(2))),
        );
        let ln = doc.get_right_neighbor(
            &get_canonical_digram(),
            &Occurrence(0, Address(AddressNumber(0), AddressNumber(1))),
        );
        assert!(ln.is_some());
        let (node_id, address_number) = ln.unwrap();
        assert_eq!(
            node_id,
            NodeId(UndirectedNodeId::new(3), Orientation::Forward)
        );
        assert_eq!(address_number, AddressNumber(2));
    }

    #[test]
    fn test_get_most_frequent() {
        let mut doc = get_doc();
        let mut occurrences: DeterministicHashSet<Occurrence> = DeterministicHashSet::default();
        occurrences.insert(Occurrence(0, Address(AddressNumber(1), AddressNumber(2))));
        occurrences.insert(Occurrence(1, Address(AddressNumber(0), AddressNumber(1))));
        doc.add_digram(&get_digram_from(2, 3).into(), occurrences);
        let freq = doc.get_most_frequent(1);
        assert!(freq.is_some());
        let freq = freq.unwrap();
        assert_eq!(freq, get_digram_from(2, 3).into());
    }
}
