use std::collections::HashMap;
use itertools::Itertools;
use std::mem;

type PathId = u64;
type Multiplicity = u32;

pub fn transpose<T>(vals: (T, T)) -> (T, T) {
    (vals.1, vals.0)
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct ColorSet(HashMap<PathId, Vec<(Multiplicity, Multiplicity)>>);

impl ColorSet {
    pub fn from(path_id: PathId, prev_count: Multiplicity, curr_count: Multiplicity) -> ColorSet {
        let mut result = ColorSet(HashMap::new());
        result.0.insert(path_id, vec![(prev_count, curr_count)]);
        result
    }

    pub fn equal_set(&self, other: &ColorSet) -> bool {
        self.0.eq(&other.0)
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn insert(&mut self, path_id: PathId, multiplicity: (Multiplicity, Multiplicity)) {
        self.0
            .entry(path_id)
            .and_modify(|v| v.push(multiplicity))
            .or_insert(vec![multiplicity]);
    }

    pub fn extend(&mut self, path_id: PathId, multiplicities: &Vec<(Multiplicity, Multiplicity)>) {
        self.0
            .entry(path_id)
            .and_modify(|v| v.extend(multiplicities.clone()))
            .or_insert(multiplicities.clone());
    }

    pub fn contains(
        &self,
        path_id: PathId,
        prev_count: Multiplicity,
        curr_count: Multiplicity,
    ) -> bool {
        self.0.get(&path_id).is_some()
            && self.0[&path_id]
                .iter()
                .any(|(a, b)| prev_count == *a && curr_count == *b)
    }

    pub fn len(&self) -> usize {
        self.0.values().flatten().count()
    }

    pub fn new() -> Self {
        Self(HashMap::new())
    }

    pub fn xu_intersection(
        &self,
        xu_set: &Self,
        mutate_outgoing: &mut HashMap<(PathId, Multiplicity), Multiplicity>,
        is_xu_flipped: bool,
        is_xq_flipped: bool,
    ) -> (ColorSet, ColorSet) {
        let mut new_xu_set = ColorSet::new();
        let mut new_xq_set = ColorSet::new();
        for path in xu_set.0.keys() {
            if let Some(uv_path_entry) = self.0.get(path) {
                for xu_step in &xu_set.0[path] {
                    let mut was_xu_step_used = false;
                    let corrected_xu_step = if is_xu_flipped {
                        transpose(*xu_step)
                    } else {
                        *xu_step
                    };
                    for uv_step in uv_path_entry {
                        if corrected_xu_step.1 == uv_step.0 {
                            was_xu_step_used = true;
                            new_xq_set.insert(
                                *path,
                                if is_xq_flipped {
                                    transpose(corrected_xu_step)
                                } else {
                                    corrected_xu_step
                                },
                            );
                            mutate_outgoing.insert((*path, uv_step.1), corrected_xu_step.1);
                        }
                    }
                    if !was_xu_step_used {
                        new_xu_set.insert(*path, *xu_step);
                    }
                }
            } else {
                new_xu_set.extend(*path, &xu_set.0[path]);
            }
        }
        (new_xu_set, new_xq_set)
    }

    pub fn vy_intersection(
        &self,
        vy_set: &Self,
        mutate_outgoing: &HashMap<(PathId, Multiplicity), Multiplicity>,
        is_vy_flipped: bool,
        is_qy_flipped: bool,
    ) -> (ColorSet, ColorSet) {
        let mut new_vy_set = ColorSet::new();
        let mut new_qy_set = ColorSet::new();
        for path in vy_set.0.keys() {
            if let Some(uv_path_entry) = self.0.get(path) {
                for vy_step in &vy_set.0[path] {
                    let mut was_vy_step_used = false;
                    let corrected_vy_step = if is_vy_flipped {
                        transpose(*vy_step)
                    } else {
                        *vy_step
                    };
                    for uv_step in uv_path_entry {
                        if corrected_vy_step.0 == uv_step.1 {
                            was_vy_step_used = true;
                            if let Some(beta) = mutate_outgoing.get(&(*path, uv_step.1)) {
                                new_qy_set.insert(
                                    *path,
                                    if is_qy_flipped {
                                        (corrected_vy_step.1, *beta)
                                    } else {
                                        (*beta, corrected_vy_step.1)
                                    },
                                )
                            } else {
                                new_qy_set.insert(
                                    *path,
                                    if is_qy_flipped {
                                        transpose(corrected_vy_step)
                                    } else {
                                        corrected_vy_step
                                    },
                                );
                            }
                        }
                    }
                    if !was_vy_step_used {
                        new_vy_set.insert(*path, *vy_step);
                    }
                }
            } else {
                new_vy_set.extend(*path, &vy_set.0[path]);
            }
        }
        (new_vy_set, new_qy_set)
    }

    pub fn self_vy_intersection(&self, mutate_outgoing: &HashMap<(PathId, Multiplicity), Multiplicity>) -> (ColorSet, ColorSet) {
        let mut new_qn_set = ColorSet::new();
        let mut new_vn_set = ColorSet::new();
        for path in self.0.keys() {
            for multiplicity in &self.0[path] {
                if let Some(replacement) = mutate_outgoing.get(&(*path, multiplicity.0)) {
                    new_qn_set.insert(*path, (*replacement, multiplicity.1));
                } else {
                    new_vn_set.insert(*path, *multiplicity);
                }
            }
        }
        (new_qn_set, new_vn_set)
    }

    pub fn sectionize(&self, mutate_outgoing: &mut HashMap<(PathId, Multiplicity), Multiplicity>) -> (ColorSet, ColorSet) {
        let (qq_set, qv_set): (HashMap<_, _>, HashMap<_, _>) = self.0.iter().map(|(path, v)| {
            let mut multiplicities = v.to_owned();
            multiplicities.sort();
            let mut sections = Vec::new();
            let mut current_section = vec![multiplicities[0]];
            multiplicities.iter().tuple_windows().for_each(|(first, second)| {
                if second.0 == first.1 || second.1 == first.0 {
                    current_section.push(*second);
                } else {
                    sections.push(mem::take(&mut current_section));
                    current_section = vec![*second];
                }
            });
            sections.push(current_section);
            let (qq_set, qv_set): (Vec<_>, Vec<_>) = sections.into_iter().map(|mut section| {
                if section[0].0 > section[0].1 {
                    section.sort_by(|a, b| b.cmp(a));
                } else {
                    section.sort();
                }
                let mut merged = Vec::new();
                let mut is_odd = false;
                section.chunks(2).for_each(|chunk| {
                    if chunk.len() == 2 {
                        merged.push((chunk[0].0, chunk[1].1));
                    } else {
                        mutate_outgoing.insert((*path, chunk[0].1), chunk[0].0);
                        is_odd = true;
                    }
                });
                let last_one = if is_odd { Vec::new() } else { merged.split_off(merged.len() - 1) };
                (merged, last_one)
            }).unzip();
            let qq_set: Vec<(u32, u32)> = qq_set.into_iter().flatten().collect();
            let qv_set: Vec<(u32, u32)> = qv_set.into_iter().flatten().collect();
            ((*path, qq_set), (*path, qv_set))
        }).unzip();

        // Remove entries for paths without multiplicities
        let qq_set = qq_set.into_iter().filter(|(_, v)| !v.is_empty()).collect();
        let qv_set = qv_set.into_iter().filter(|(_, v)| !v.is_empty()).collect();

        (ColorSet(qq_set), ColorSet(qv_set))
    }
}

impl PartialOrd for ColorSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(
            self.0
                .iter()
                .map(|(_, v)| v)
                .flatten()
                .count()
                .cmp(&other.0.iter().map(|(_, v)| v).flatten().count()),
        )
    }
}

impl Ord for ColorSet {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.len().cmp(&other.0.len())
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_self_vy_contained() {
        let c = ColorSet(HashMap::from([(0, vec![(3,3)])]));
        let mutate_outgoing = HashMap::from([((0, 3), 2)]);
        let (qn_set, vn_set) = c.self_vy_intersection(&mutate_outgoing);
        assert_eq!(qn_set, ColorSet(HashMap::from([(0, vec![(2,3)])])));
        assert_eq!(vn_set, ColorSet(HashMap::new()));
    }

    #[test]
    fn test_self_vy_not_contained() {
        let c = ColorSet(HashMap::from([(0, vec![(3,3)])]));
        let mutate_outgoing = HashMap::new();
        let (qn_set, vn_set) = c.self_vy_intersection(&mutate_outgoing);
        assert_eq!(qn_set, ColorSet(HashMap::new()));
        assert_eq!(vn_set, ColorSet(HashMap::from([(0, vec![(3,3)])])));
    }

    #[test]
    fn test_sectionize_inc_odd() {
        let c = ColorSet(HashMap::from([(0, vec![(1,2), (2,3), (3,6), (6,7), (7,8)])]));
        let mut mutate_outgoing = HashMap::new();
        let (qq_set, qv_set) = c.sectionize(&mut mutate_outgoing);
        assert_eq!(qq_set, ColorSet(HashMap::from([(0, vec![(1,3), (3,7)])])));
        assert_eq!(qv_set, ColorSet::new());
        assert_eq!(mutate_outgoing, HashMap::from([((0, 8), 7)]));
    }

    #[test]
    fn test_sectionize_inc_even() {
        let c = ColorSet(HashMap::from([(0, vec![(1,2), (2,3), (3,6), (6,7), (7,8), (8,9)])]));
        let mut mutate_outgoing = HashMap::new();
        let (qq_set, qv_set) = c.sectionize(&mut mutate_outgoing);
        assert_eq!(qq_set, ColorSet(HashMap::from([(0, vec![(1,3), (3,7)])])));
        assert_eq!(qv_set, ColorSet(HashMap::from([(0, vec![(7,9)])])));
        assert_eq!(mutate_outgoing, HashMap::new());
    }

    #[test]
    fn test_sectionize_dec_odd() {
        let c = ColorSet(HashMap::from([(0, vec![(8,7), (7,6), (6,3), (3,2), (2,1)])]));
        let mut mutate_outgoing = HashMap::new();
        let (qq_set, qv_set) = c.sectionize(&mut mutate_outgoing);
        assert_eq!(qq_set, ColorSet(HashMap::from([(0, vec![(8,6), (6,2)])])));
        assert_eq!(qv_set, ColorSet::new());
        assert_eq!(mutate_outgoing, HashMap::from([((0, 1), 2)]));
    }

    #[test]
    fn test_sectionize_dec_even() {
        let c = ColorSet(HashMap::from([(0, vec![(9,8), (8,7), (7,6), (6,3), (3,2), (2,1)])]));
        let mut mutate_outgoing = HashMap::new();
        let (qq_set, qv_set) = c.sectionize(&mut mutate_outgoing);
        assert_eq!(qq_set, ColorSet(HashMap::from([(0, vec![(9,7), (7,3)])])));
        assert_eq!(qv_set, ColorSet(HashMap::from([(0, vec![(3,1)])])));
        assert_eq!(mutate_outgoing, HashMap::new());
    }

    #[test]
    fn test_sectionize_multiple() {
        let c = ColorSet(HashMap::from([(0, vec![(9,8), (8,7), (7,6), (1,2), (2,3)])]));
        let mut mutate_outgoing = HashMap::new();
        let (qq_set, qv_set) = c.sectionize(&mut mutate_outgoing);
        assert_eq!(qq_set, ColorSet(HashMap::from([(0, vec![(9,7)])])));
        assert_eq!(qv_set, ColorSet(HashMap::from([(0, vec![(1,3)])])));
        assert_eq!(mutate_outgoing, HashMap::from([((0, 6), 7)]));

    }

    #[test]
    fn test_xu_1() {
        let mut xu = ColorSet::from(0, 1, 2);
        xu.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let mut mutate_outgoing = HashMap::new();
        let (new_xu, new_xq) = uv.xu_intersection(&xu, &mut mutate_outgoing, false, false);
        assert_eq!(new_xu, ColorSet::from(0, 2, 3));
        assert_eq!(new_xq, ColorSet::from(0, 1, 2));
        assert_eq!(mutate_outgoing, HashMap::from([((0, 3), 2)]));
    }

    #[test]
    fn test_xu_2() {
        let mut xu = ColorSet::from(0, 2, 1);
        xu.insert(0, (3, 2));
        let uv = ColorSet::from(0, 2, 3);
        let mut mutate_outgoing = HashMap::new();
        let (new_xu, new_xq) = uv.xu_intersection(&xu, &mut mutate_outgoing, true, false);
        assert_eq!(new_xu, ColorSet::from(0, 3, 2));
        assert_eq!(new_xq, ColorSet::from(0, 1, 2));
    }

    #[test]
    fn test_xu_3() {
        let mut xu = ColorSet::from(0, 1, 2);
        xu.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let mut mutate_outgoing = HashMap::new();
        let (new_xu, new_xq) = uv.xu_intersection(&xu, &mut mutate_outgoing, false, true);
        assert_eq!(new_xu, ColorSet::from(0, 2, 3));
        assert_eq!(new_xq, ColorSet::from(0, 2, 1));
    }

    #[test]
    fn test_xu_4() {
        let mut xu = ColorSet::from(0, 2, 1);
        xu.insert(0, (3, 2));
        let uv = ColorSet::from(0, 2, 3);
        let mut mutate_outgoing = HashMap::new();
        let (new_xu, new_xq) = uv.xu_intersection(&xu, &mut mutate_outgoing, true, true);
        assert_eq!(new_xu, ColorSet::from(0, 3, 2));
        assert_eq!(new_xq, ColorSet::from(0, 2, 1));
    }

    #[test]
    fn test_vy_1() {
        let mut vy = ColorSet::from(0, 3, 4);
        vy.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &HashMap::new(), false, false);
        assert_eq!(new_vy, ColorSet::from(0, 2, 3));
        assert_eq!(new_qy, ColorSet::from(0, 3, 4));
    }

    #[test]
    fn test_vy_2() {
        let mut vy = ColorSet::from(0, 4, 3);
        vy.insert(0, (3, 2));
        let uv = ColorSet::from(0, 2, 3);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &HashMap::new(), true, false);
        assert_eq!(new_vy, ColorSet::from(0, 3, 2));
        assert_eq!(new_qy, ColorSet::from(0, 3, 4));
    }

    #[test]
    fn test_vy_3() {
        let mut vy = ColorSet::from(0, 3, 4);
        vy.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &HashMap::new(), false, true);
        assert_eq!(new_vy, ColorSet::from(0, 2, 3));
        assert_eq!(new_qy, ColorSet::from(0, 4, 3));
    }

    #[test]
    fn test_vy_4() {
        let mut vy = ColorSet::from(0, 4, 3);
        vy.insert(0, (3, 2));
        let uv = ColorSet::from(0, 2, 3);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &HashMap::new(), true, true);
        assert_eq!(new_vy, ColorSet::from(0, 3, 2));
        assert_eq!(new_qy, ColorSet::from(0, 4, 3));
    }

    #[test]
    fn test_vy_substitution_1() {
        let mut vy = ColorSet::from(0, 3, 4);
        vy.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let mutate_outgoing = HashMap::from([((0, 3), 7)]);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &mutate_outgoing, false, false);
        assert_eq!(new_vy, ColorSet::from(0, 2, 3));
        assert_eq!(new_qy, ColorSet::from(0, 7, 4));
    }

    #[test]
    fn test_vy_substitution_2() {
        let mut vy = ColorSet::from(0, 3, 4);
        vy.insert(0, (2, 3));
        let uv = ColorSet::from(0, 2, 3);
        let mutate_outgoing = HashMap::from([((0, 3), 7)]);
        let (new_vy, new_qy) = uv.vy_intersection(&vy, &mutate_outgoing, false, true);
        assert_eq!(new_vy, ColorSet::from(0, 2, 3));
        assert_eq!(new_qy, ColorSet::from(0, 4, 7));
    }
}
