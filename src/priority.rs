use crate::node_id::NodeId;
use std::collections::{HashMap, HashSet};

pub struct Priority {
    rules: HashMap<NodeId, FRule>,
    priority_list: Vec<HashSet<NodeId>>,
    index_lookup: HashMap<i64, usize>,
}

impl Priority {
    pub fn into_rules(self) -> HashMap<NodeId, Vec<NodeId>> {
        self.rules.into_iter().map(|(k, v)| v.into_rule()).collect()
    }

    pub fn from(rules: HashMap<NodeId, FRule>) -> Self {
        let mut scores: HashSet<i64> = HashSet::new();
        for (rule_name, rule) in &rules {
            scores.insert(rule.get_score());
        }
        let mut scores: Vec<i64> = scores.into_iter().collect();
        scores.sort();
        let mut priority_list: Vec<HashSet<NodeId>> = vec![HashSet::new(); scores.len()];
        let index_lookup: HashMap<i64, usize> = scores.iter().enumerate().map(|(i, s)| (*s, i)).collect();
        for (rule_name, rule) in &rules {
            let score = rule.get_score();
            let idx = index_lookup[&score];
            priority_list[idx].insert(*rule_name);
        }
        Self {
            rules,
            priority_list,
            index_lookup,
        }
    }

    pub fn get_max_rule(&self) -> Option<NodeId> {
        let mut counter = self.priority_list.len() - 1;
        loop {
            if !self.priority_list[counter].is_empty() {
                let first = self.priority_list[counter].iter().next().expect("Has at least one element");
                return Some(first.to_owned());
            }
            counter -= 1;
            if counter == 0 {
                return None;
            }
        }
    }

    fn find_closest(&self, new_score: i64) -> (i64, usize) {
        let mut closest_score = 0i64;
        let mut closest_idx = 0usize;
        for (score, idx) in &self.index_lookup {
            if *score > closest_score && *score  < new_score {
                closest_score = *score;
                closest_idx = *idx;
            }
        }
        (closest_score, closest_idx)
    }

    fn insert_new_score(&mut self, new_score: i64, closest_idx: usize) {
        self.index_lookup.iter_mut().for_each(|(k, v)| if *v >= (closest_idx + 1) { *v += 1 });
        self.index_lookup.insert(new_score, closest_idx + 1);
    }

    fn second_branch(&mut self, key: &NodeId, new_score: i64) {
        let (closest_score, closest_idx) = self.find_closest(new_score);
        let new_idx = closest_idx + 1;
        self.insert_new_score(new_score, closest_idx);
        self.priority_list.insert(closest_idx + 1, HashSet::from([*key]));
    }

    pub fn change_occurrence(&mut self, key: &NodeId, occurrence: usize) {
        let old_score = self.rules[key].get_score();
        self.rules.get_mut(key).expect("Priority rules has key").occurrence = occurrence;
        let new_score = self.rules[key].get_score();
        if new_score < 0 {
            panic!("Negative score for {}: {} {} -> {}", key, occurrence, old_score, new_score);
        }
        let old_idx = self.index_lookup[&old_score];
        self.priority_list[old_idx].remove(key);

        if self.priority_list[old_idx].is_empty() && (self.priority_list.len() - 1) == old_idx {
            self.index_lookup.remove(&old_score);
            self.priority_list.remove(old_idx);
        }
        let branch = self.index_lookup.contains_key(&new_score);
        let mut new_idx = 0;
        if self.index_lookup.contains_key(&new_score) {
            new_idx = self.index_lookup[&new_score];
            self.priority_list[new_idx].insert(*key);
        } else {
            // Find closest without going over
            self.second_branch(key, new_score);
        }
        //let mut wrong_ones = vec![];
        //for (rule_name, rule) in &self.rules {
        //    let score = rule.get_score();
        //    let idx = self.index_lookup[&score];
        //    if self.priority_list[idx].iter().position(|x| *x == *rule_name).is_none() {
        //        wrong_ones.push((*rule_name, idx));
        //    }
        //}
        //if !wrong_ones.is_empty() {
        //    eprintln!("Something went wrong during changeOccurrence {}: {}, branch: {}, new_idx: {}", key, occurrence, branch, new_idx);
        //    eprintln!("not at correct position:");
        //    for (el, s) in &wrong_ones {
        //        eprintln!("\t{}: {}", el, s);
        //    }
        //    panic!();
        //}
    }

    pub fn subtract_occurrence(&mut self, key: &NodeId, value: usize) {
        let old_occurrance = self.rules[key].occurrence;
        self.change_occurrence(key, old_occurrance - value);
    }

    pub fn get_children(&self, key: &NodeId) -> Vec<NodeId> {
        self.rules[key].right.clone()
    }

    pub fn get_parents(&self, key: &NodeId) -> Vec<NodeId> {
        self.rules[key].parents.clone()
    }

    pub fn get_full_rule(&self, rule: &NodeId) -> Vec<NodeId> {
        let children = self.get_children(rule);
        let mut full_rule: Vec<NodeId> = Vec::new();
        for oriented_child in &children {
            let child = &oriented_child.get_forward();
            if self.rules.contains_key(child) {
                let full_sub_rule = self.get_full_rule(child);
                let full_sub_rule = if !oriented_child.is_forward() { full_sub_rule.into_iter().rev().map(|x| x.flip()).collect::<Vec<NodeId>>() } else { full_sub_rule };
                full_rule.extend(full_sub_rule);
            } else {
                full_rule.push(*oriented_child);
            }
        }
        full_rule
    }

    pub fn len(&self) -> usize {
        self.rules.len()
    }

    pub fn get_occurrence(&self, key: &NodeId) -> usize {
        self.rules[key].occurrence
    }

    pub fn contains_key(&self, key: &NodeId) -> bool {
        self.rules.contains_key(key)
    }
}

#[derive(Clone, Debug)]
pub struct FRule {
    pub left: NodeId,
    pub right: Vec<NodeId>,
    pub length: usize,
    pub occurrence: usize,
    pub parents: Vec<NodeId>,
}

impl FRule {
    pub fn into_rule(self) -> (NodeId, Vec<NodeId>) {
        (self.left, self.right)
    }

    pub fn get_score(&self) -> i64 {
        let score = (self.length as i64 - 1) * self.occurrence as i64 - self.length as i64 - 1;
        if score < 0 {
            0
        } else {
            score
        }
    }
}
