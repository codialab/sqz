use std::collections::{HashMap, VecDeque};
use std::fmt::Debug;
use std::hash::Hash;

use itertools::Itertools;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match<T> {
    pub pattern_id: T,
    pub pattern_range: std::ops::Range<usize>,
}

struct TrieNode<T> {
    next_node: HashMap<T, usize>,
    fail: usize,
    matches: Vec<(T, usize)>,
}

impl<T> TrieNode<T> {
    fn new() -> Self {
        Self {
            next_node: HashMap::new(),
            fail: 0,
            matches: Vec::new(),
        }
    }
}

pub struct AhoCorasick<T> {
    nodes: Vec<TrieNode<T>>,
}

impl<T: Hash + Eq + Clone + Debug> AhoCorasick<T> {
    pub fn num_states(&self) -> usize {
        self.nodes.len()
    }

    pub fn new(patterns: &HashMap<T, Vec<T>>) -> Self {
        let mut ac = Self {
            nodes: vec![TrieNode::new()],
        };

        for (k, v) in patterns.iter() {
            ac.insert(v, k);
        }

        ac.build_failure_links();

        ac
    }

    fn insert(&mut self, pattern: &[T], id: &T) {
        let mut curr = 0;
        let pattern_length = pattern.len();
        let iter = pattern.iter();

        for token in iter {
            if self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].next_node[token];
            } else {
                self.nodes.push(TrieNode::new());
                let last_node = self.nodes.len() - 1;
                self.nodes
                    .get_mut(curr)
                    .unwrap()
                    .next_node
                    .insert(token.clone(), last_node);
                curr = last_node;
            }
        }
        self.nodes
            .get_mut(curr)
            .unwrap()
            .matches
            .push((id.clone(), pattern_length));
    }

    fn build_failure_links(&mut self) {
        let mut queue: VecDeque<usize> = VecDeque::new();
        for (_, child) in self.nodes[0].next_node.iter() {
            queue.push_back(*child);
        }

        while let Some(curr_parent) = queue.pop_front() {
            for (token, child) in self.nodes[curr_parent].next_node.clone() {
                queue.push_back(child);
                let mut parent_of_fail = self.nodes[curr_parent].fail;
                while parent_of_fail != 0
                    && !self.nodes[parent_of_fail].next_node.contains_key(&token)
                {
                    parent_of_fail = self.nodes[parent_of_fail].fail;
                }
                let failure_target = match self.nodes[parent_of_fail].next_node.get(&token) {
                    Some(fail_node) => *fail_node,
                    None => 0,
                };
                self.nodes[child].fail = failure_target;
                let matches_of_failure_target = self.nodes[failure_target].matches.clone();
                self.nodes[child].matches.extend(matches_of_failure_target);
            }
        }
    }

    #[allow(dead_code)]
    pub fn find_all(&self, haystack: &[T]) -> Vec<Match<T>> {
        let mut matches = Vec::new();
        let mut curr = 0;

        let mut cache = Vec::new();

        for (pos, token) in haystack.iter().enumerate() {
            while curr != 0 && !self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].fail;
            }
            if self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].next_node[token];
                cache.extend(
                    self.nodes[curr]
                        .matches
                        .iter()
                        .map(|(pattern, length)| Match {
                            pattern_id: pattern.clone(),
                            pattern_range: (pos + 1 - length)..(pos + 1),
                        }),
                );
            }
        }

        if cache.is_empty() {
            return matches;
        }

        cache.sort_by(|a, b| {
            a.pattern_range
                .start
                .cmp(&b.pattern_range.start)
                .then(b.pattern_range.len().cmp(&a.pattern_range.len()))
        });

        let mut current_end = 0;
        for candidate in cache {
            if candidate.pattern_range.start >= current_end {
                current_end = candidate.pattern_range.end;
                matches.push(candidate);
            }
        }

        matches
    }
}

impl<T: Hash + Eq + Clone + Debug + ReverseComplementable> AhoCorasick<T> {
    pub fn find_all_reverse_complement(&self, haystack: &[T]) -> Vec<Match<T>> {
        let mut matches = Vec::new();
        let mut curr = 0;

        let mut cache = Vec::new();

        for (pos, token) in haystack.iter().enumerate() {
            while curr != 0 && !self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].fail;
            }
            if self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].next_node[token];
                cache.extend(
                    self.nodes[curr]
                        .matches
                        .iter()
                        .map(|(pattern, length)| Match {
                            pattern_id: pattern.clone(),
                            pattern_range: (pos + 1 - length)..(pos + 1),
                        }),
                );
            }
        }

        let reverse_haystack = haystack.iter().rev().map(|t| t.complement()).collect_vec();
        let length_haystack = reverse_haystack.len();
        let mut curr = 0;
        for (rev_pos, token) in reverse_haystack.iter().enumerate() {
            while curr != 0 && !self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].fail;
            }
            if self.nodes[curr].next_node.contains_key(token) {
                curr = self.nodes[curr].next_node[token];
                cache.extend(
                    self.nodes[curr]
                        .matches
                        .iter()
                        .map(|(pattern, length)| Match {
                            pattern_id: pattern.complement(),
                            pattern_range: (length_haystack - (rev_pos + 1))
                                ..(length_haystack - (rev_pos + 1) + length),
                        }),
                );
            }
        }

        if cache.is_empty() {
            return matches;
        }

        cache.sort_by(|a, b| {
            a.pattern_range
                .start
                .cmp(&b.pattern_range.start)
                .then(b.pattern_range.len().cmp(&a.pattern_range.len()))
        });

        let mut current_end = 0;
        for candidate in cache {
            if candidate.pattern_range.start >= current_end {
                current_end = candidate.pattern_range.end;
                matches.push(candidate);
            }
        }

        matches
    }
}

pub trait ReverseComplementable {
    fn complement(&self) -> Self;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::helpers::utils::simple_node_creation::*;

    #[test]
    fn test_simple_search() {
        let patterns = HashMap::from([(0, vec![1, 2]), (1, vec![2, 3, 4])]);
        let ac = AhoCorasick::new(&patterns);
        // Haystack: 0, 1, 2, 3, 4, 5
        let haystack = vec![1, 2, 3, 4];

        let matches = ac.find_all(&haystack);

        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0],
            Match {
                pattern_id: 0,
                pattern_range: 0..2
            }
        )
    }

    #[test]
    fn test_find_longest() {
        let patterns = HashMap::from([(0, vec![1, 2]), (1, vec![1, 2, 3])]);
        let ac = AhoCorasick::new(&patterns);
        // Haystack: 0, 1, 2, 3, 4, 5
        let haystack = vec![1, 2, 3, 4];

        let matches = ac.find_all(&haystack);

        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0],
            Match {
                pattern_id: 1,
                pattern_range: 0..3
            }
        )
    }

    #[test]
    fn test_find_rev_comp() {
        let patterns = HashMap::from([
            (m(10), vec![n(1), n(2), n(3)]),
            (m(11), vec![rn(4), rn(3), rn(2), rn(1)]),
        ]);
        let ac = AhoCorasick::new(&patterns);
        let haystack = vec![n(1), n(2), n(3), n(4)];
        let matches = ac.find_all_reverse_complement(&haystack);

        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0],
            Match {
                pattern_id: rm(11),
                pattern_range: 0..4
            }
        )
    }
}
