use itertools::Itertools;

use crate::{
    grammar_building::Rule,
    helpers::{
        DeterministicHashMap, digram_occurrences::DigramOccurrences, utils::{Address, Digram, NodeId}
    },
};

pub fn get_haplotype_walks(
    d: &DigramOccurrences,
    rules: &[Rule],
    singleton_haplotypes: &DeterministicHashMap<usize, NodeId>,
    number_of_paths: usize,
) -> Vec<Vec<NodeId>> {
    let mut walks: Vec<Vec<(Digram, Address)>> = vec![Vec::new(); number_of_paths];

    // Find all singleton digrams
    for (digram, occurrences) in d {
        for occurrence in occurrences {
            let walk_id = occurrence.get_haplotype();
            let address = occurrence.get_address();
            if address.is_forward() {
                walks[walk_id].push((digram.clone().into(), address));
            } else {
                walks[walk_id].push((Into::<Digram>::into(digram.clone()).flip(), address.flip()));
            }
        }
    }

    log::debug!("walks: {:?}", walks);

    // Sort all singleton digrams
    let mut walks = walks
        .into_iter()
        .map(|mut walk| {
            walk.sort_by(|a, b| a.1.cmp(&b.1));
            if !walk.is_empty() {
                let first_value = vec![walk[0].0 .0];
                first_value
                    .into_iter()
                    .chain(walk.into_iter().map(|(digram, _address)| digram.1))
                    .collect_vec()
            } else {
                Vec::new()
            }
        })
        .collect_vec();

    // Insert all walks that only consist of a single node
    let single_node_walks = walks
        .iter()
        .enumerate()
        .filter_map(|(idx, walk)| if walk.is_empty() { Some(idx) } else { None })
        .collect_vec();
    log::info!("SNW: {:?}", single_node_walks);
    for single_node_walk in single_node_walks {
        if singleton_haplotypes.contains_key(&single_node_walk) {
            walks[single_node_walk].push(singleton_haplotypes[&single_node_walk]);
        } else if let Some(meta_node) = get_single_meta_node_walk(single_node_walk, rules) {
            walks[single_node_walk].push(meta_node);
        } else {
            log::error!(
                "Was not able to assign a sequence to walk No. {}",
                single_node_walk
            );
        }
    }
    walks
}

fn get_single_meta_node_walk(single_node_walk: usize, rules: &[Rule]) -> Option<NodeId> {
    for rule in rules.iter().rev() {
        for occurrence in &rule.3 {
            if occurrence.get_haplotype() == single_node_walk {
                if occurrence.get_address().is_forward() {
                    return Some(rule.0);
                } else {
                    return Some(rule.0.flip());
                }
            }
        }
    }
    None
}
