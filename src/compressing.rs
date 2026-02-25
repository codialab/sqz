use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;

use itertools::Itertools;
use rayon::prelude::*;

use crate::aho_corasick::AhoCorasick;
use crate::aho_corasick::Match;
use crate::helpers::PathSegment;
use crate::{
    helpers::{utils::NodeId, DeterministicHashMap, NodeRegistry, ReverseNodeRegistry},
    parser::{bufreader_from_compressed, path_from_line, ByteLineReader},
};

pub fn compress_remaining_file(
    file: &PathBuf,
    mut rules: DeterministicHashMap<NodeId, Vec<NodeId>>,
    compressed_paths: &[bool],
    node_reg: &NodeRegistry,
    rev_reg: &ReverseNodeRegistry,
) -> Result<()> {
    unroll_rules(&mut rules);
    log::info!("Unrolled all rules");

    log::info!("Collected all small rules");

    let ac = {
        let patterns: HashMap<NodeId, Vec<NodeId>> = rules.into_iter().collect();

        AhoCorasick::new(&patterns)
    };
    log::info!(
        "Done with setting up Aho-Corasick data structure with {} states",
        ac.num_states()
    );

    let data = bufreader_from_compressed(file)?;
    let line_reader = ByteLineReader::new(data);
    line_reader
        .filter_map(|line| path_from_line(&line, node_reg))
        .enumerate()
        .filter(|(idx, _)| !compressed_paths[*idx])
        .batching(|it| {
            let mut batch = Vec::with_capacity(50);
            for _ in 0..50 {
                match it.next() {
                    Some(item) => batch.push(item),
                    None => break,
                }
            }
            if batch.is_empty() {
                None
            } else {
                Some(batch)
            }
        })
        .par_bridge()
        .for_each(|c| {
            for (_, (name, sequence)) in c {
                let compressed_text = compress_haplotype(sequence, &ac);
                let compressed_text: String = compressed_text
                    .into_iter()
                    .map(|n| rev_reg.get_directed_name(n))
                    .collect();
                println!("W\t{}\t{}", name, compressed_text);
            }
        });
    log::info!("Done with compressing remaining file");
    Ok(())
}

pub fn compress(
    haplotypes: Vec<(PathSegment, Vec<NodeId>)>,
    mut rules: DeterministicHashMap<NodeId, Vec<NodeId>>,
) -> Vec<(PathSegment, Vec<NodeId>)> {
    unroll_rules(&mut rules);
    log::info!("Unrolled all rules");
    let ac = {
        let patterns: HashMap<NodeId, Vec<NodeId>> = rules.into_iter().collect();

        AhoCorasick::new(&patterns)
    };
    log::info!(
        "Done with setting up Aho-Corasick data structure with {} states",
        ac.num_states()
    );
    haplotypes
        .into_iter()
        .map(|(n, h)| (n, compress_haplotype(h, &ac)))
        .collect()
}

// Unrolls rules so that no rule text contains meta-nodes
fn unroll_rules(rules: &mut DeterministicHashMap<NodeId, Vec<NodeId>>) {
    let keys: Vec<NodeId> = rules.keys().copied().collect();
    for meta_node in keys {
        unroll_rule(&meta_node, rules);
    }
}

fn unroll_rule(key: &NodeId, rules: &mut DeterministicHashMap<NodeId, Vec<NodeId>>) {
    if rules[key].iter().all(|n| !n.is_meta_node()) {
        return;
    }

    let mut new_seq: Vec<NodeId> = Vec::new();
    let old_seq = rules[key].clone();
    for node in old_seq {
        if node.is_meta_node() {
            unroll_rule(&node.get_forward(), rules);
            let seq = if node.is_forward() {
                rules[&node].clone()
            } else {
                rules[&node.get_forward()]
                    .iter()
                    .rev()
                    .map(|n| n.flip())
                    .collect()
            };
            new_seq.extend(seq.into_iter());
        } else {
            new_seq.push(node);
        }
    }
    *rules.get_mut(key).expect("Can find key meta-node") = new_seq;
}

fn compress_haplotype(haplotype: Vec<NodeId>, ac: &AhoCorasick<NodeId>) -> Vec<NodeId> {
    let matches = ac.find_all_reverse_complement(&haplotype);
    let mut text = haplotype.into_iter().map(Some).collect_vec();
    for Match {
        pattern_id,
        pattern_range,
    } in matches
    {
        let start = pattern_range.start;
        let end = pattern_range.end;
        text[start] = Some(pattern_id);
        for entry in text.iter_mut().take(end).skip(start + 1) {
            *entry = None;
        }
    }

    text.into_iter().flatten().collect_vec()
}

#[cfg(test)]
mod tests {
    use crate::helpers::utils::UndirectedNodeId;
    use crate::parser::parse_walk_seq_plainly;

    use super::*;

    fn get_node(id: u32, meta: bool, orientation: bool) -> NodeId {
        let undir = UndirectedNodeId(id, meta);
        NodeId(
            undir,
            if orientation {
                crate::helpers::utils::Orientation::Forward
            } else {
                crate::helpers::utils::Orientation::Backward
            },
        )
    }

    #[test]
    fn test_comress_simple() {
        let h = parse_walk_seq_plainly(">1>2>3>4>5".as_bytes());
        let mut rules = DeterministicHashMap::default();
        rules.insert(
            get_node(6, true, true),
            vec![get_node(1, false, true), get_node(2, false, true)],
        );
        rules.insert(
            get_node(7, true, true),
            vec![get_node(3, false, true), get_node(4, false, true)],
        );
        rules.insert(
            get_node(8, true, true),
            vec![
                get_node(6, true, true),
                get_node(7, true, true),
                get_node(5, false, true),
            ],
        );
        let mut res = compress(
            vec![(
                PathSegment::new(
                    "A".to_string(),
                    "0".to_string(),
                    "1".to_string(),
                    None,
                    None,
                ),
                h,
            )],
            rules,
        );
        let res = std::mem::take(&mut res[0].1);
        assert_eq!(res, vec![get_node(8, true, true)]);
    }
}
